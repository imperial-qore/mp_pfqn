#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <gmp.h>
#include <gmpla.h>
#include "util.h"
#include "safe_tiers.h"

/* Shared implementation of the tiered "always exact" wrappers; see
 * safe_tiers.h for the tier structure and safe_comom/README.md for the
 * relation to TSE'09 Appendix A. */

static char BINDIR[4096];

/* run "<BINDIR>/<solver> <model> <flags>"; capture stdout into out (cap
 * bytes); set *singular if the run reports a singular system or fails.
 * Returns number of bytes read, or -1 on spawn failure. */
static int run_solver(const char* wrapper, const char* solver, const char* model,
                      const char* flags, char* out, int cap, int* singular)
{
	char cmd[8192], errf[256];
	/* stdout carries the clean numeric result; stderr carries progress
	 * and the singular/perturbation messages.  Keep them separate: parse
	 * stdout for output, both streams for the singular check. */
	snprintf(errf, sizeof(errf), "/tmp/%s_err_%d", wrapper, (int)getpid());
	snprintf(cmd, sizeof(cmd), "%s/%s '%s' %s 2>%s ; echo EXIT:$?",
	         BINDIR, solver, model, flags ? flags : "", errf);
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

/* "<flag> <solver_flags>", the invocation used for Tier 1 and for the
 * final Tier-2 solve. */
static void join_flags(char* dst, int cap, const char* flag, const char* extra)
{
	if (extra && extra[0]) snprintf(dst, cap, "%s %s", flag ? flag : "", extra);
	else                   snprintf(dst, cap, "%s", flag ? flag : "");
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
 * -t has R rows of 1 value (one throughput per line); -q has M rows of R
 * values. */
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
 * class -- O(R) orders, not O(R!) -- probing singularity with the cheap
 * probe flags; solve and print the first non-singular one.  Returns 1 on
 * success. */
static int recclass_search(const safe_cfg* cfg, qnmodel* qn, const char* flag,
                           int is_t, int is_q)
{
	/* a wrapper whose outputs carry no class index needs no un-permutation */
	const char* pf = cfg->permuted_flags ? cfg->permuted_flags : "-t -q";
	if (!strstr(pf, "-t")) is_t = 0;
	if (!strstr(pf, "-q")) is_q = 0;
	int R = qn->R, c, j;
	int perm[64];
	char run_flags[256];
	join_flags(run_flags, sizeof(run_flags), flag, cfg->solver_flags);
	for (c = 0; c < R; c++) {
		if (c == R-1) continue;               /* identity already tried (Tier 1) */
		/* class c last, the others in original order */
		int p = 0;
		for (j = 0; j < R; j++) if (j != c) perm[p++] = j;
		perm[R-1] = c;
		char tmp[128]; snprintf(tmp, sizeof(tmp), "/tmp/%s_%d.qn", cfg->name, (int)getpid());
		write_perm(tmp, qn, perm);
		char out[1<<20]; int sing;
		run_solver(cfg->name, cfg->solver, tmp, cfg->probe_flags, out, sizeof(out), &sing);
		if (sing) { unlink(tmp); continue; }
		run_solver(cfg->name, cfg->solver, tmp, run_flags, out, sizeof(out), &sing);
		unlink(tmp);
		if (!sing) {
			if (is_t || is_q) print_unpermuted(out, R, qn->M, perm, is_q);
			else fputs(out, stdout);
			return 1;
		}
	}
	return 0;
}

int safe_tiers_main(const safe_cfg* cfg, int argc, char** argv)
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
		printf("USAGE: %s %s model.qn\n", argv[0],
		       cfg->flags_usage ? cfg->flags_usage : "[-e|-g|-l|-t|-q]");
		if (cfg->usage) fputs(cfg->usage, stdout);
		return -1;
	}

	/* locate the sibling binaries (bin/ dir of argv[0]) */
	char self[4096]; strncpy(self, argv[0], sizeof(self)-1); self[sizeof(self)-1] = '\0';
	strncpy(BINDIR, dirname(self), sizeof(BINDIR)-1);

	char out[1<<20]; int sing;
	char run_flags[256];
	join_flags(run_flags, sizeof(run_flags), flag, cfg->solver_flags);

	/* Tier 1: the primary solver as given */
	run_solver(cfg->name, cfg->solver, model, run_flags, out, sizeof(out), &sing);
	if (!sing) { fputs(out, stdout); return 0; }

	/* Tier 2: the primary solver over class permutations */
	qnmodel* qn = readmodel((char*)model);
	if (qn->hasOpen) {
		/* open/mixed models are outside the closed-network domain of the
		 * Tier-3 convolution; report rather than emit a bogus result. */
		fputs(out, stdout);
		fprintf(stderr, "%s: open/mixed model is out of scope\n", cfg->name);
		return 1;
	}
	if (!qn->isLD) {
		if (recclass_search(cfg, qn, flag, is_t, is_q)) return 0;
	}

	if (!cfg->fallback) {
		/* No exact Tier 3 exists for this output (see safe_tiers.h).  Say
		 * which of the two reasons applies: a load-dependent model never
		 * reached the class search at all, so calling it singular would be
		 * wrong. */
		if (qn->isLD)
			fprintf(stderr, "%s: load-dependent model is outside %s's domain\n",
			        cfg->name, cfg->solver);
		else
			fprintf(stderr, "%s: %s singular under every class order and no exact\n"
			                "  fallback exists for this output; no answer\n",
			        cfg->name, cfg->solver);
		return 1;
	}

	/* Tier 3: exact convolution (always exact) */
	fprintf(stderr, "%s: %s singular under every class order; using exact %s (Tier 3)\n",
	        cfg->name, cfg->solver, cfg->fallback);
	run_solver(cfg->name, cfg->fallback, model, flag, out, sizeof(out), &sing);
	if (sing) { fprintf(stderr, "%s: %s fallback failed\n", cfg->name, cfg->fallback); return 1; }
	fputs(out, stdout);
	return 0;
}
