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

/* Output modes that carry per-class values (ported from main's safe_comom).
 * OUT_RAW: class-order invariant, printed verbatim. */
enum { OUT_RAW = 0, OUT_T, OUT_Q, OUT_TEXT };

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
				printmodel(qn);          /* same routine the solvers use, original order */
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
 * class -- O(R) orders, not O(R!) -- probing singularity with the cheap
 * probe flags; solve and print the first non-singular one.  Returns 1 on
 * success. */
static int recclass_search(const safe_cfg* cfg, qnmodel* qn, const char* flag,
                           int is_t, int is_q)
{
	/* a wrapper whose outputs carry no class index needs no un-permutation;
	 * "-d" in permuted_flags also remaps the default/-d text report */
	const char* pf = cfg->permuted_flags ? cfg->permuted_flags : "-t -q";
	int is_text = (!flag[0] || !strcmp(flag, "-d") || !strcmp(flag, "--debug"));
	int mode = OUT_RAW;
	if (is_t && strstr(pf, "-t")) mode = OUT_T;
	else if (is_q && strstr(pf, "-q")) mode = OUT_Q;
	else if (is_text && strstr(pf, "-d")) mode = OUT_TEXT;
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
			if (mode == OUT_RAW) fputs(out, stdout);
			else if (print_unpermuted(out, qn, perm, mode) < 0) {
				fprintf(stderr, "%s: unexpected %s output under class order; trying the next one\n",
				        cfg->name, cfg->solver);
				continue;
			}
			return 1;
		}
	}
	return 0;
}

int safe_tiers_main(const safe_cfg* cfg, int argc, char** argv)
{
	const char* flag = ""; int is_t = 0, is_q = 0, nflags = 0;
	char* model = NULL; int i;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			/* the solvers print only one output, so several flags would be
			 * misread when un-permuting a Tier-2 run */
			if (++nflags > 1) {
				fprintf(stderr, "%s: give at most one output option (got '%s' and '%s')\n",
				        cfg->name, flag, argv[i]);
				return 1;
			}
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
