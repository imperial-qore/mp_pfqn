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
 * -t (per class), -q (per station,class), -d and the default report from a
 * permuted Tier-2 run are mapped back to the original class order before
 * printing.  Open and load-dependent models are refused up front.
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

/* Output modes of comom that carry per-class values. */
enum { OUT_G = 0, OUT_T, OUT_Q, OUT_TEXT };

/* print a result captured from a permuted run, mapped back to the original
 * class order.  perm[j] = original class at permuted position j.
 *   OUT_T    : R lines, one throughput each
 *   OUT_Q    : M lines of R values
 *   OUT_TEXT : default / -d report.  The model echo (printed in permuted
 *              order) is replaced by printmodel() of the original model, the
 *              "X[j] = v" lines and the values of the "Q[k] =\tv..\t(total)"
 *              lines are reordered; everything else passes through.
 * Returns 0 on success, -1 if the output does not have the expected shape
 * (nothing is printed in that case). */
static int print_unpermuted(const char* out, qnmodel* qn, int* perm, int mode)
{
	int R = qn->R, M = qn->M, i;
	int* inv = (int*)malloc((size_t)R * sizeof(int)); /* inv[orig] = permuted position */
	for (i = 0; i < R; i++) inv[perm[i]] = i;

	if (mode == OUT_T || mode == OUT_Q) {
		int cap = (mode == OUT_Q ? M*R : R), nv = 0;
		char** vals = (char**)calloc((size_t)cap, sizeof(char*));
		char* buf = strdup(out);
		char* tok = strtok(buf, " \t\n");
		while (tok && nv < cap) { vals[nv++] = tok; tok = strtok(NULL, " \t\n"); }
		if (nv != cap || tok != NULL) { free(vals); free(buf); free(inv); return -1; }
		if (mode == OUT_Q) {
			int k, c;
			for (k = 0; k < M; k++) {
				for (c = 0; c < R; c++)
					printf("%s%s", vals[k*R + inv[c]], c<R-1?" ":"");
				printf("\n");
			}
		} else {
			int c;
			for (c = 0; c < R; c++) printf("%s\n", vals[inv[c]]);
		}
		free(vals); free(buf); free(inv);
		return 0;
	}

	/* OUT_TEXT: split into lines, rewrite into a buffer, print at the end */
	char* buf = strdup(out);
	int nlines = 0, cap_lines = 64;
	char** line = (char**)malloc((size_t)cap_lines * sizeof(char*));
	{
		char* p = buf;
		while (*p) {
			char* e = strchr(p, '\n');
			if (nlines == cap_lines) { cap_lines *= 2; line = (char**)realloc(line, (size_t)cap_lines * sizeof(char*)); }
			line[nlines++] = p;
			if (!e) break;
			*e = '\0'; p = e + 1;
		}
	}
	/* X lines: the R consecutive "X[j] = " lines; values in permuted order */
	int x0 = -1;
	for (i = 0; i < nlines; i++) if (!strncmp(line[i], "X[", 2)) { x0 = i; break; }
	int ok = (x0 >= 0 && x0 + R <= nlines);
	for (i = 0; ok && i < R; i++) {
		char want[32]; snprintf(want, sizeof(want), "X[%d] = ", i+1);
		if (strncmp(line[x0+i], want, strlen(want))) ok = 0;
	}
	/* Q lines: "Q[k] =" then R tab-separated values then "\t(total: ..)" */
	int q0 = -1;
	for (i = 0; i < nlines; i++) if (!strncmp(line[i], "Q[", 2)) { q0 = i; break; }
	if (ok && !(q0 >= 0 && q0 + M <= nlines)) ok = 0;
	char*** qv = NULL;
	if (ok) {
		qv = (char***)calloc((size_t)M, sizeof(char**));
		for (i = 0; ok && i < M; i++) {
			int nt = 0;
			qv[i] = (char**)calloc((size_t)R + 2, sizeof(char*));
			char* t = strtok(line[q0+i], "\t");
			while (t && nt < R + 2) { qv[i][nt++] = t; t = strtok(NULL, "\t"); }
			if (nt != R + 2 || t != NULL || strncmp(qv[i][R+1], "(total:", 7)) ok = 0;
		}
	}
	if (ok) {
		int in_header = 0;
		for (i = 0; i < nlines; i++) {
			if (!in_header && !strncmp(line[i], "The queueing network model has", 30)) {
				in_header = 1;
				fflush(stdout);
				printmodel(qn);          /* same routine comom uses, original order */
				fflush(stdout);
				continue;
			}
			if (in_header) {             /* the echo ends at its blank line */
				if (line[i][0] == '\0') { in_header = 0; printf("\n"); }
				continue;
			}
			if (i >= x0 && i < x0 + R) {
				int c = i - x0;          /* original class c, permuted position inv[c] */
				printf("X[%d] = %s\n", c+1, line[x0 + inv[c]] + strcspn(line[x0 + inv[c]], "=") + 2);
			} else if (i >= q0 && i < q0 + M) {
				int k = i - q0, c;
				printf("%s", qv[k][0]);
				for (c = 0; c < R; c++) printf("\t%s", qv[k][1 + inv[c]]);
				printf("\t%s\n", qv[k][R+1]);
			} else {
				int last_nl = (i < nlines-1) || (out[0] && out[strlen(out)-1] == '\n');
				printf("%s%s", line[i], last_nl ? "\n" : "");
			}
		}
	}
	if (qv) { for (i = 0; i < M; i++) free(qv[i]); free(qv); }
	free(line); free(buf); free(inv);
	return ok ? 0 : -1;
}

/* Recursion-class search (the appendix's reorder idea).  The recursion
 * class's loadings do not enter the coefficient matrix, so a singularity
 * that depends on which class is processed last is removed by moving a
 * suitable class to that position.  Try each class as the last (recursion)
 * class -- O(R) orders, not O(R!) -- probing singularity with the fast
 * comom -e path; solve and print the first non-singular one.  Returns 1 on
 * success. */
static int recclass_search(qnmodel* qn, const char* flag, int mode)
{
	int R = qn->R, c, j;
	int* perm = (int*)malloc((size_t)R * sizeof(int));
	static char out[1<<20];
	for (c = 0; c < R; c++) {
		if (c == R-1) continue;               /* identity already tried (Tier 1) */
		/* class c last, the others in original order */
		int p = 0;
		for (j = 0; j < R; j++) if (j != c) perm[p++] = j;
		perm[R-1] = c;
		char tmp[64]; snprintf(tmp, sizeof(tmp), "/tmp/safe_comom_%d.qn", (int)getpid());
		write_perm(tmp, qn, perm);
		int sing;
		run_solver("comom", tmp, "-e", out, sizeof(out), &sing);   /* fast probe */
		if (sing) { unlink(tmp); continue; }
		run_solver("comom", tmp, flag, out, sizeof(out), &sing);
		unlink(tmp);
		if (!sing) {
			if (mode == OUT_G) fputs(out, stdout);
			else if (print_unpermuted(out, qn, perm, mode) < 0) {
				fprintf(stderr, "safe_comom: unexpected comom output under class order; trying the next one\n");
				continue;
			}
			free(perm);
			return 1;
		}
	}
	free(perm);
	return 0;
}

int main(int argc, char** argv)
{
	/* At most one output option: comom prints only one of them, so passing
	 * several would make the class un-permutation below misread its output. */
	const char* flag = ""; int mode = OUT_TEXT, nflags = 0;
	char* model = NULL; int i;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			const char* a = argv[i];
			int m;
			if (!strcmp(a,"-e") || !strcmp(a,"--ex") || !strcmp(a,"-g") || !strcmp(a,"--nc") ||
			    !strcmp(a,"-l") || !strcmp(a,"--log")) m = OUT_G;
			else if (!strcmp(a,"-t") || !strcmp(a,"--tput")) m = OUT_T;
			else if (!strcmp(a,"-q") || !strcmp(a,"--qlen")) m = OUT_Q;
			else if (!strcmp(a,"-d") || !strcmp(a,"--debug")) m = OUT_TEXT;
			else {
				fprintf(stderr, "safe_comom: unsupported option '%s' (use one of -e -g -l -t -q -d)\n", a);
				return 1;
			}
			if (++nflags > 1) {
				fprintf(stderr, "safe_comom: give at most one output option (got '%s' and '%s')\n", flag, a);
				return 1;
			}
			flag = a; mode = m;
		} else model = argv[i];
	}
	if (!model) {
		printf("USAGE: %s [-e|-g|-l|-t|-q|-d] model.qn\n", argv[0]);
		printf("  CoMoM that always returns the exact answer, via the hybrid\n");
		printf("  MVA/CoMoM of TSE'09 Appendix A: plain CoMoM, then CoMoM with a\n");
		printf("  permuted class order, then exact convolution for genuine\n");
		printf("  loading degeneracies.\n");
		return -1;
	}

	/* locate the sibling binaries (bin/ dir of argv[0]) */
	char self[4096]; strncpy(self, argv[0], sizeof(self)-1); self[sizeof(self)-1] = '\0';
	strncpy(BINDIR, dirname(self), sizeof(BINDIR)-1);

	static char out[1<<20]; int sing;

	/* Open classes and load-dependent (MU) stations are outside CoMoM's
	 * domain.  Refuse them up front: comom used to ignore MU and print a
	 * wrong answer that Tier 1 accepted, and a refusal from comom must not
	 * be mistaken for a singularity and sent on to convolution. */
	qnmodel* qn = readmodel((char*)model);
	qn_require_closed(qn, "safe_comom", 0);

	/* Tier 1: plain CoMoM */
	run_solver("comom", model, flag, out, sizeof(out), &sing);
	if (!sing) { fputs(out, stdout); return 0; }

	/* Tier 2: CoMoM over class permutations */
	if (recclass_search(qn, flag, mode)) return 0;

	/* Tier 3: exact convolution (always exact) */
	fprintf(stderr, "safe_comom: CoMoM singular under every class order; using exact convolution (Tier 3)\n");
	run_solver("ca", model, flag, out, sizeof(out), &sing);
	if (sing) { fprintf(stderr, "safe_comom: convolution fallback failed\n"); return 1; }
	fputs(out, stdout);
	return 0;
}
