#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <gmp.h>
#include <gmpla.h>
#include "util.h"

/* safe_comom - CoMoM that always returns an exact answer, using the
 * hybrid MVA/CoMoM idea of Casale, "CoMoM" (IEEE TSE 2009), Appendix A:
 * when the CoMoM coefficient matrix is singular for a degenerate model,
 * interleave convolution with CoMoM so the degenerate directions are
 * handled by an MVA/convolution-style recursion.
 *
 * Realised as a tiered orchestrator over the exact solvers already in the
 * tree; every tier is exact, so the result is always the exact G/X/Q.
 *
 *   Tier 1  plain CoMoM (bin/comom).
 *   Tier 2  CoMoM after permuting the class order.  A singularity that
 *           depends only on which class is the recursion class (the last
 *           one) is removed by processing a non-degenerate class last;
 *           this is the cheap case of the appendix (permute class indices).
 *   Tier 3  exact convolution (bin/ca) for genuine loading degeneracies
 *           that no class order removes.  This is the fully-peeled limit
 *           of the appendix's MVA-like recursion (all degenerate classes
 *           swept by convolution); it is always exact.
 *
 * G(N) is invariant to the class order, so -e/-g/-l need no un-permutation;
 * -t (per class) and -q (per station,class) from a permuted Tier-2 run are
 * mapped back to the original class order before printing.
 *
 * The appendix's efficiency refinement -- peeling only the MINIMAL set of
 * degenerate classes by convolution while keeping CoMoM for the rest,
 * rather than falling all the way back to convolution in Tier 3 -- is a
 * further optimisation over this exact-but-tiered version; see README.
 */

static char BINDIR[4096];

/* run "<BINDIR>/<solver> <model> <flag>"; capture stdout into out (cap
 * bytes); set *singular if the run reports a singular system or fails.
 * Returns number of bytes read, or -1 on spawn failure. */
static int run_solver(const char* solver, const char* model, const char* flag,
                      char* out, int cap, int* singular)
{
	char cmd[8192], errf[128];
	/* stdout carries the clean numeric result; stderr carries progress
	 * and the singular/perturbation messages.  Keep them separate: parse
	 * stdout for output, both streams for the singular check. */
	snprintf(errf, sizeof(errf), "/tmp/safe_comom_err_%d", (int)getpid());
	snprintf(cmd, sizeof(cmd), "%s/%s '%s' %s 2>%s ; echo EXIT:$?",
	         BINDIR, solver, model, flag ? flag : "", errf);
	FILE* p = popen(cmd, "r");
	if (!p) return -1;
	int n = 0, c;
	while ((c = fgetc(p)) != EOF && n < cap-1) out[n++] = (char)c;
	out[n] = '\0';
	pclose(p);
	int rc = 0; char* e = strstr(out, "EXIT:");
	if (e) { rc = atoi(e+5); *e = '\0'; n = (int)(e - out); }

	/* read captured stderr */
	char err[8192]; err[0] = '\0';
	FILE* ef = fopen(errf, "r");
	if (ef) { int m = fread(err, 1, sizeof(err)-1, ef); err[m>0?m:0] = '\0'; fclose(ef); }
	unlink(errf);

	/* a run is unusable if it failed or reported a singular system (the
	 * exact solvers either print an error, or warn and fall back to
	 * perturbation); reject both in favour of the next tier. */
	#define HITS(b) (strstr((b),"singular")||strstr((b),"Singular")|| \
	                 strstr((b),"perturbation")||strstr((b),"Perturbation")|| \
	                 strstr((b),"Error")||strstr((b),"NaN")||strstr((b),"nan"))
	*singular = (rc != 0) || HITS(out) || HITS(err);
	#undef HITS
	return n;
}

/* write a class-permuted copy of the model to path.  perm[j] = original
 * class placed at permuted position j. */
static void write_perm(const char* path, qnmodel* qn, int* perm)
{
	FILE* f = fopen(path, "w");
	int R = qn->R, M = qn->M, j, k;
	fprintf(f, "%d\n", R);
	for (j = 0; j < R; j++) gmp_fprintf(f, "%d%s", qn->N[perm[j]], j<R-1?" ":"\n");
	for (j = 0; j < R; j++) gmp_fprintf(f, "%Zd%s", qn->Z[perm[j]], j<R-1?" ":"\n");
	fprintf(f, "%d\n", M);
	for (k = 0; k < M; k++) {
		fprintf(f, "%d", qn->mi[k]);
		for (j = 0; j < R; j++) gmp_fprintf(f, " %Zd", qn->L[k][perm[j]]);
		fprintf(f, "\n");
	}
	fclose(f);
}

/* print a -t/-q result captured from a permuted run, mapped back to the
 * original class order.  perm[j]=orig class at permuted position j.
 * ncols_per_row = R for -t rows (1 row) ... actually -t has R rows of 1
 * value (one throughput per line); -q has M rows of R values. */
static void print_unpermuted(const char* out, int R, int M, int* perm, int is_q)
{
	/* inverse permutation: inv[orig] = permuted position */
	int inv[64]; int i;
	for (i = 0; i < R; i++) inv[perm[i]] = i;

	/* parse all whitespace-separated doubles */
	int cap = (is_q ? M*R : R);
	char* vals[4096]; int nv = 0;
	char* buf = strdup(out);
	char* tok = strtok(buf, " \t\n");
	while (tok && nv < cap) { vals[nv++] = strdup(tok); tok = strtok(NULL, " \t\n"); }
	free(buf);

	if (is_q) {
		int k, c;
		for (k = 0; k < M; k++) {
			for (c = 0; c < R; c++) {
				printf("%s%s", vals[k*R + inv[c]], c<R-1?" ":"");
			}
			printf("\n");
		}
	} else {
		int c;
		for (c = 0; c < R; c++) printf("%s\n", vals[inv[c]]);
	}
	for (i = 0; i < nv; i++) free(vals[i]);
}

/* Recursion-class search (the appendix's reorder idea).  The recursion
 * class's loadings do not enter the coefficient matrix, so a singularity
 * that depends on which class is processed last is removed by moving a
 * suitable class to that position.  Try each class as the last (recursion)
 * class -- O(R) orders, not O(R!) -- probing singularity with the fast
 * comom -e path; solve and print the first non-singular one.  Returns 1 on
 * success. */
static int recclass_search(qnmodel* qn, const char* flag, int is_t, int is_q)
{
	int R = qn->R, c, j;
	int perm[64];
	for (c = 0; c < R; c++) {
		if (c == R-1) continue;               /* identity already tried (Tier 1) */
		/* class c last, the others in original order */
		int p = 0;
		for (j = 0; j < R; j++) if (j != c) perm[p++] = j;
		perm[R-1] = c;
		char tmp[64]; snprintf(tmp, sizeof(tmp), "/tmp/safe_comom_%d.qn", (int)getpid());
		write_perm(tmp, qn, perm);
		char out[1<<20]; int sing;
		run_solver("comom", tmp, "-e", out, sizeof(out), &sing);   /* fast probe */
		if (sing) { unlink(tmp); continue; }
		run_solver("comom", tmp, flag, out, sizeof(out), &sing);
		unlink(tmp);
		if (!sing) {
			if (is_t || is_q) print_unpermuted(out, R, qn->M, perm, is_q);
			else fputs(out, stdout);
			return 1;
		}
	}
	return 0;
}

int main(int argc, char** argv)
{
	const char* flag = ""; int is_t = 0, is_q = 0;
	char* model = NULL; int i;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			flag = argv[i];
			if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tput")) is_t = 1;
			if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) is_q = 1;
		} else model = argv[i];
	}
	if (!model) {
		printf("USAGE: %s [-e|-g|-l|-t|-q] model.qn\n", argv[0]);
		printf("  CoMoM that always returns the exact answer, via the hybrid\n");
		printf("  MVA/CoMoM of TSE'09 Appendix A: plain CoMoM, then CoMoM with a\n");
		printf("  permuted class order, then exact convolution for genuine\n");
		printf("  loading degeneracies.\n");
		return -1;
	}

	/* locate the sibling binaries (bin/ dir of argv[0]) */
	char self[4096]; strncpy(self, argv[0], sizeof(self)-1); self[sizeof(self)-1] = '\0';
	strncpy(BINDIR, dirname(self), sizeof(BINDIR)-1);

	char out[1<<20]; int sing;

	/* Tier 1: plain CoMoM */
	run_solver("comom", model, flag, out, sizeof(out), &sing);
	if (!sing) { fputs(out, stdout); return 0; }

	/* Tier 2: CoMoM over class permutations */
	qnmodel* qn = readmodel((char*)model);
	if (qn->hasOpen) {
		/* open/mixed models are outside CoMoM's (and ca's) closed-network
		 * domain; report rather than emit a bogus convolution result. */
		fputs(out, stdout);
		fprintf(stderr, "safe_comom: open/mixed model is out of scope for CoMoM\n");
		return 1;
	}
	if (!qn->isLD) {
		if (recclass_search(qn, flag, is_t, is_q)) return 0;
	}

	/* Tier 3: exact convolution (always exact) */
	fprintf(stderr, "safe_comom: CoMoM singular under every class order; using exact convolution (Tier 3)\n");
	run_solver("ca", model, flag, out, sizeof(out), &sing);
	if (sing) { fprintf(stderr, "safe_comom: convolution fallback failed\n"); return 1; }
	fputs(out, stdout);
	return 0;
}
