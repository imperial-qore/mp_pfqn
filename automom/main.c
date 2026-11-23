/* automom - pick and run the fastest exact solver of the MoM family for a model.
 *
 * For each candidate solver the work is estimated in closed form from its
 * own basis sizes (exact formulas, checked against binmom/bincomom -v):
 *   mom basis   d_m(r) = r [C(M+r-1,r) + C(M+r-2,r-1)]
 *   comom basis d_c(r) = (M+1) C(M+r-1,r-1)
 *   gcomom      U_r = sum over its sub-model family of |support| C(D+r-1,r-1)
 *   gmom        colA(m,r) = r C(m+r-2,r-1) at prefix m (normal equations every step)
 *   ca          prod_r (N_r+1) R sum_k mi_k
 * and the bit size b of G (about Ntot log2(Lmax sum mi + Zmax)) for the cost of
 * one rational operation.  The predicted time is exp(c0) W^c1 b^c2 with
 * (c0, c1, c2) fitted per solver by automom/calibrate.py (calib.h).
 *
 * The candidates allowed by the model (domain, flags) are ranked by predicted
 * time and run in that order until one succeeds: a refusal (binmom on mixed
 * think times or a singular system, gmom -X on a singular system) or a
 * failure moves on to the next.  Load-dependent and mixed models, outside the
 * MoM family, go to comomld/gld or mvaldmx/mvamx. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <gmp.h>
#include "util.h"
#include "calib.h"

enum { S_CA, S_MOM, S_COMOM, S_GCOMOM, S_GMOM, S_BINMOM, S_BINCOMOM, S_MOMMOD, S_COMOMMOD, NSOLV };
static const char* NAME[NSOLV] = { "ca", "mom", "comom", "gcomom", "gmom", "binmom", "bincomom", "mommod", "comommod" };
static const bool HAS_J[NSOLV] = { false, true, true, false, false, true, true, true, true };

static double lnck(int n, int k) { if (k < 0 || n < k) return -INFINITY; return lgamma(n+1.0) - lgamma(k+1.0) - lgamma(n-k+1.0); }
static double nckd(int n, int k) { return exp(lnck(n, k)); }

/* ---- gcomom family: copy vectors closed under removing one copy of each of
 * the top min(R-1, |support|) support stations (same rule as gcomom) ---- */
typedef struct { int M, n, cap, hcap; int* v; int* key; } fam;
static uint64_t code(const int* v, int M) { uint64_t h = 1469598103934665603ULL; for (int k = 0; k < M; k++) { h ^= (uint64_t)(v[k]+1); h *= 1099511628211ULL; } return h; }
static int fam_add(fam* F, const int* v)
{
	uint64_t h = code(v, F->M) & (F->hcap-1);
	while (F->key[h] >= 0) { if (!memcmp(F->v + (size_t)F->key[h]*F->M, v, F->M*sizeof(int))) return 0; h = (h+1) & (F->hcap-1); }
	if (F->n == F->cap) { F->cap *= 2; F->v = realloc(F->v, (size_t)F->cap*F->M*sizeof(int)); }
	memcpy(F->v + (size_t)F->n*F->M, v, F->M*sizeof(int));
	F->key[h] = F->n++;
	if (2*F->n > F->hcap) {
		free(F->key); F->hcap *= 2; F->key = malloc(F->hcap*sizeof(int));
		for (int z = 0; z < F->hcap; z++) F->key[z] = -1;
		for (int z = 0; z < F->n; z++) { uint64_t g = code(F->v + (size_t)z*F->M, F->M) & (F->hcap-1); while (F->key[g] >= 0) g = (g+1) & (F->hcap-1); F->key[g] = z; }
	}
	return 1;
}

/* U[r] for r = 1..R; returns the family size, or -1 above cap */
static long gcomom_unknowns(qnmodel* q, double* U, long cap)
{
	int M = q->M, R = q->R;
	fam F = { M, 0, 1024, 4096, NULL, NULL };
	F.v = malloc((size_t)F.cap*M*sizeof(int)); F.key = malloc(F.hcap*sizeof(int));
	for (int z = 0; z < F.hcap; z++) F.key[z] = -1;
	int* t = malloc(M*sizeof(int));
	fam_add(&F, q->mi);
	for (int r = 1; r <= R; r++) U[r] = 0;
	long n = 0;
	for (int id = 0; id < F.n; id++) {
		if (F.n > cap) { n = -1; break; }
		int* v = F.v + (size_t)id*M, nt = 0;
		for (int k = 0; k < M; k++) if (v[k] > 0) nt++;
		int D = nt >= R-1 ? 1 : (nt > 0 ? nt : 1);
		for (int r = 1; r <= R; r++) U[r] += nt * nckd(D + r - 1, r - 1);
		int nb = MIN(R-1, nt), seen = 0;
		for (int k = M-1; k >= 0 && seen < nb; k--) {
			if (v[k] <= 0) continue;
			seen++;
			memcpy(t, F.v + (size_t)id*M, M*sizeof(int)); t[k]--;
			int any = 0; for (int z = 0; z < M; z++) if (t[z]) { any = 1; break; }
			if (any) fam_add(&F, t);
			v = F.v + (size_t)id*M;   /* fam_add may move F.v */
		}
	}
	if (n == 0) n = F.n;
	free(F.v); free(F.key); free(t);
	return n;
}

/* ---- model features ---- */
typedef struct {
	int M, R, Ntot; double logb;
	bool ld, open, mi1, zmixed, distinct_rows;
	double W[NSOLV];          /* structural work, INFINITY = not estimated */
	long famsize;
} feat;

static void features(qnmodel* q, feat* f)
{
	int M = q->M, R = q->R;
	memset(f, 0, sizeof *f);
	f->M = M; f->R = R;
	f->ld = q->isLD; f->open = q->hasOpen;
	for (int r = 0; r < R; r++) { if (q->N[r] < 0) f->open = true; else f->Ntot += q->N[r]; }
	f->mi1 = true; long smi = 0;
	for (int k = 0; k < M; k++) { if (q->mi[k] != 1) f->mi1 = false; smi += q->mi[k]; }
	int nz = 0, zz = 0; for (int r = 0; r < R; r++) { if (mpz_sgn(q->Z[r])) nz++; else zz++; }
	f->zmixed = nz > 0 && zz > 0;
	double Lmax = 1, Zmax = 0;
	for (int k = 0; k < M; k++) for (int r = 0; r < R; r++) Lmax = fmax(Lmax, mpz_get_d(q->L[k][r]));
	for (int r = 0; r < R; r++) Zmax = fmax(Zmax, mpz_get_d(q->Z[r]));
	f->logb = log(f->Ntot * log2(2.0 + Lmax*smi + Zmax) + 64.0);
	if (f->ld || f->open) { for (int s = 0; s < NSOLV; s++) f->W[s] = INFINITY; return; }

	double ca = R * (double)smi; for (int r = 0; r < R; r++) ca *= (q->N[r] + 1.0);
	double mom = 0, binmom = 0, comom = 0, bincomom = 0, gmom = 0;
	for (int r = 1; r <= R; r++) {
		int Nr = q->N[r-1];
		double dm = r * (nckd(M+r-1, r) + nckd(M+r-2, r-1));
		double dc = (M + 1) * nckd(M+r-1, r-1);
		double sm = dm*dm*M, ss = (double)Nr * dm * (M+1);
		mom += sm + ss;
		binmom += sm + ((Nr >= 2*dm*dm) ? dm*dm*dm*log2(Nr+1.0) : ss);
		comom += dc*dc*dc + (double)Nr*dc*dc;
		bincomom += dc*dc*dc + ((Nr - M - 1 >= 2*dc*dc) ? dc*dc*dc*log2(Nr+1.0) : (double)Nr*dc*dc);
		if (r >= 2) for (int m = 2; m <= M; m++) { double ca_ = r * nckd(m+r-2, r-1); gmom += (Nr + 1.0) * ca_*ca_*ca_; }
	}
	double U[64]; double gc = 0;
	f->famsize = (R < 63) ? gcomom_unknowns(q, U, 400000) : -1;
	if (f->famsize > 0) for (int r = 1; r <= R; r++) gc += U[r] * (q->N[r-1] + (double)R*R) * R;
	double np = 2.0 * exp(f->logb) / 61.0 + 2;    /* primes for num and den */
	f->W[S_CA] = ca; f->W[S_MOM] = mom; f->W[S_COMOM] = comom;
	f->W[S_GCOMOM] = f->famsize > 0 ? gc : INFINITY;
	f->W[S_GMOM] = gmom > 0 ? gmom : INFINITY;
	f->W[S_BINMOM] = binmom; f->W[S_BINCOMOM] = bincomom;
	f->W[S_MOMMOD] = mom * np; f->W[S_COMOMMOD] = comom * np;
}

static bool allowed(const feat* f, int s, bool dbg)
{
	if (!isfinite(f->W[s])) return false;
	switch (s) {
	case S_GMOM:   return f->R >= 2 && f->mi1 && !dbg;
	case S_BINMOM: return !f->zmixed;
	case S_MOMMOD: return !dbg;          /* no -d */
	default:       return true;
	}
}

static double predict(const feat* f, int s, int threads)
{
	double t = exp(CALIB[s].c0 + CALIB[s].c1 * log(fmax(f->W[s], 1.0)) + CALIB[s].c2 * f->logb);
	if (threads > 1 && HAS_J[s] && (s == S_MOMMOD || s == S_COMOMMOD)) t /= threads;
	return t + 0.003;
}

/* ---- running a solver ---- */
static char BINDIR[PATH_MAX];

static int run(const char* solver, char** flags, int nflags, const char* jarg, const char* model, bool extra_X)
{
	char path[PATH_MAX + 64];
	snprintf(path, sizeof path, "%s/%s", BINDIR, solver);
	char* av[32]; int n = 0;
	av[n++] = (char*) solver;
	for (int i = 0; i < nflags; i++) av[n++] = flags[i];
	if (extra_X) av[n++] = "-X";
	if (jarg) { av[n++] = "-j"; av[n++] = (char*) jarg; }
	av[n++] = (char*) model; av[n] = NULL;
	fflush(stdout); fflush(stderr);
	pid_t p = fork();
	if (p == 0) { execv(path, av); fprintf(stderr, "automom: cannot run %s\n", path); _exit(127); }
	int st; waitpid(p, &st, 0);
	return WIFEXITED(st) ? WEXITSTATUS(st) : 128;
}

/* buffered run: the output of a failed attempt is discarded */
static int run_capture(const char* solver, char** flags, int nflags, const char* jarg, const char* model, bool extra_X, bool quiet)
{
	char tmpl[] = "/tmp/automom.XXXXXX";
	int fd = mkstemp(tmpl);
	if (fd < 0) return run(solver, flags, nflags, jarg, model, extra_X);
	fflush(stdout);
	int saved = dup(1); dup2(fd, 1);
	int saved_err = -1, devnull = -1;
	if (quiet) { saved_err = dup(2); devnull = open("/dev/null", 1); dup2(devnull, 2); }
	int rc = run(solver, flags, nflags, jarg, model, extra_X);
	fflush(stdout); dup2(saved, 1); close(saved);
	if (quiet) { dup2(saved_err, 2); close(saved_err); close(devnull); }
	if (rc == 0) {
		lseek(fd, 0, SEEK_SET);
		char buf[65536]; ssize_t k;
		while ((k = read(fd, buf, sizeof buf)) > 0) fwrite(buf, 1, k, stdout);
	}
	close(fd); unlink(tmpl);
	return rc;
}

static void usage(const char* p)
{
	printf("Usage: %s [options] model.qn\n", p);
	printf("Runs the exact solver of the MoM family predicted to be fastest for the model.\n");
	printf("  -l -e -g -t -q   output, as in every solver (-d is passed through; its format is the solver's)\n");
	printf("  -j n             threads, for the solvers that take -j\n");
	printf("  --explain        print the candidates, their predicted times and the choice on stderr\n");
	printf("  --dry-run        print the choice and exit\n");
	printf("  --solver NAME    run NAME (one of ca mom comom gcomom gmom binmom bincomom mommod comommod)\n");
	printf("  --features       print the cost features (used by automom/calibrate.py)\n");
}

int main(int argc, char** argv)
{
	char* flags[16]; int nf = 0;
	bool explain = false, dry = false, featmode = false, dbg = false, nc_only = true;
	const char* jarg = NULL; const char* force = NULL; char* file = NULL;
	for (int i = 1; i < argc; i++) {
		const char* a = argv[i];
		if (!strcmp(a, "-l") || !strcmp(a, "-e") || !strcmp(a, "-g") || !strcmp(a, "-t") || !strcmp(a, "-q") || !strcmp(a, "-d")) {
			if (nf < 15) flags[nf++] = argv[i];
			if (!strcmp(a, "-d")) dbg = true;
			if (!strcmp(a, "-t") || !strcmp(a, "-q") || !strcmp(a, "-d")) nc_only = false;
		}
		else if (!strcmp(a, "-j") && i+1 < argc) jarg = argv[++i];
		else if (!strcmp(a, "--explain")) explain = true;
		else if (!strcmp(a, "--dry-run")) { dry = true; explain = true; }
		else if (!strcmp(a, "--features")) featmode = true;
		else if (!strcmp(a, "--solver") && i+1 < argc) force = argv[++i];
		else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0]); return 0; }
		else if (a[0] == '-') { fprintf(stderr, "automom: unknown option %s\n", a); return 1; }
		else file = argv[i];
	}
	if (!file) { usage(argv[0]); return 1; }
	if (nf == 0) nc_only = false;
	/* solvers live next to this binary */
	char self[PATH_MAX]; ssize_t L = readlink("/proc/self/exe", self, sizeof self - 1);
	if (L > 0) { self[L] = 0; snprintf(BINDIR, sizeof BINDIR, "%s", dirname(self)); } else snprintf(BINDIR, sizeof BINDIR, "./bin");

	qnmodel* q = readmodel(file);
	if (!q) { fprintf(stderr, "automom: cannot read %s\n", file); return 1; }
	nckinit(MAXNCKTABLE - 1, MAXNCKTABLE - 1);
	feat f; features(q, &f);
	int threads = jarg ? atoi(jarg) : 1;

	if (featmode) {
		printf("M %d R %d Ntot %d logb %.6f famsize %ld\n", f.M, f.R, f.Ntot, f.logb, f.famsize);
		for (int s = 0; s < NSOLV; s++) printf("%s %.6e %d\n", NAME[s], f.W[s], allowed(&f, s, false));
		return 0;
	}

	/* outside the MoM family */
	if (f.ld || f.open) {
		const char* seq[3] = { NULL, NULL, NULL };
		if (f.open) seq[0] = f.ld ? "mvaldmx" : "mvamx";
		else if (nc_only) { seq[0] = "comomld"; seq[1] = "gld"; }
		else seq[0] = "mvaldmx";
		if (explain) fprintf(stderr, "automom: %s%s model, outside the MoM family: %s%s%s\n",
		                     f.open ? "mixed" : "closed", f.ld ? " load-dependent" : "", seq[0], seq[1] ? ", then " : "", seq[1] ? seq[1] : "");
		if (dry) return 0;
		for (int i = 0; i < 2 && seq[i]; i++) {
			int rc = (seq[i+1]) ? run_capture(seq[i], flags, nf, NULL, file, false, false) : run(seq[i], flags, nf, NULL, file, false);
			if (rc == 0) return 0;
		}
		return 1;
	}

	/* rank */
	int ord[NSOLV], n = 0; double pt[NSOLV];
	for (int s = 0; s < NSOLV; s++) {
		pt[s] = allowed(&f, s, dbg) ? predict(&f, s, threads) : INFINITY;
		if (force) { if (!strcmp(force, NAME[s])) ord[n++] = s; }
		else if (isfinite(pt[s])) ord[n++] = s;
	}
	if (force && n == 0) { fprintf(stderr, "automom: unknown solver %s\n", force); return 1; }
	for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) if (pt[ord[j]] < pt[ord[i]]) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
	if (explain) {
		fprintf(stderr, "automom: M=%d R=%d Ntot=%d, gcomom family %ld sub-models\n", f.M, f.R, f.Ntot, f.famsize);
		for (int s = 0; s < NSOLV; s++) {
			if (isfinite(pt[s])) fprintf(stderr, "  %-9s predicted %10.3g s\n", NAME[s], pt[s]);
			else fprintf(stderr, "  %-9s %s\n", NAME[s], isfinite(f.W[s]) ? "not applicable (domain or flags)" : "not estimated");
		}
		fprintf(stderr, "automom: running %s\n", NAME[ord[0]]);
	}
	if (dry) return 0;
	for (int i = 0; i < n; i++) {
		int s = ord[i];
		const char* j = (jarg && HAS_J[s]) ? jarg : NULL;
		bool X = (s == S_GMOM);            /* never let gmom fall back to an approximate perturbation */
		int rc = (i < n-1) ? run_capture(NAME[s], flags, nf, j, file, X, s == S_GMOM || s == S_BINMOM)
		                   : run(NAME[s], flags, nf, j, file, X);
		if (rc == 0) return 0;
		if (i < n-1) fprintf(stderr, "automom: %s failed or declined (exit %d), trying %s\n", NAME[s], rc, NAME[ord[i+1]]);
	}
	return 1;
}
