/* automom - pick and run the fastest exact solver of the MoM family for a model.
 *
 * For each candidate solver the work is estimated in closed form from its
 * own basis sizes (exact formulas, checked against binmom/bincomom -v):
 *   mom basis   d_m(r) = r [C(M+r-1,r) + C(M+r-2,r-1)]
 *   comom basis d_c(r) = (M+1) C(M+r-1,r-1)
 *   gcomom      U_r = sum over its sub-model family of |support| C(D+r-1,r-1)
 *   gmom        sum over its model tree (Bmax = 1 for gmom, R-1 for gmomb) of
 *               colA^3 + N_r colA^2, colA = r C(d+l-1,l), l = max(1, r-min(d,Bmax))
 *   ca          prod_r (N_r+1) R sum_k mi_k
 * and the bit size b of G for the cost of one rational operation, from the
 * a-priori bounds of util/gbound.h (log2 G_up plus the denominator bound),
 * which also give the exact prime count of mommod/comommod.  The predicted time is exp(c0) W^c1 b^c2 with
 * (c0, c1, c2) fitted per solver by automom/calibrate.py (calib.h).
 *
 * The candidates allowed by the model (domain, flags) are ranked by predicted
 * time and run in that order until one succeeds: a refusal (binmom on mixed
 * think times or a singular system, gmom -X on a singular system) or a
 * failure moves on to the next.  --race runs the first two at once when
 * their predictions are close and keeps the first to succeed.  Load-dependent and mixed models, outside the
 * MoM family, go to comomld/gld or mvaldmx/mvamx; a fractional multiplicity
 * (util.h qn_mult) only to mom, comom, gcomom and gmom/gmomb. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <sys/wait.h>
#include <signal.h>
#include <fcntl.h>
#include <gmp.h>
#include "util.h"
#include "gbound.h"
#include "calib.h"
#include "calib_frac.h"   /* constants for fractional multiplicities (calibrate.py --frac) */

enum { S_CA, S_MOM, S_COMOM, S_GCOMOM, S_GMOM, S_BINMOM, S_BINCOMOM, S_MOMMOD, S_COMOMMOD, S_GMOMB, NSOLV };
static const char* NAME[NSOLV] = { "ca", "mom", "comom", "gcomom", "gmom", "binmom", "bincomom", "mommod", "comommod", "gmomb" };
/* gmomb is bin/gmom --bmax R-1, the branching factor that brings every
 * sub-model with R-1 or more queues to basis level 1 (PEVA'11) */
static const char* BIN[NSOLV]  = { "ca", "mom", "comom", "gcomom", "gmom", "binmom", "bincomom", "mommod", "comommod", "gmom" };
static const bool HAS_J[NSOLV] = { false, true, true, true, true, true, true, true, true, true };

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

/* fractional multiplicity mi/mden (util.h): a never-removed base fb[k] =
 * mi % mden plus vfull[k] = mi / mden removable copies, as in gcomom and
 * gmom; integer models have fb = 0, vfull = mi */
static void split_mi(qnmodel* q, int* fb, int* vfull)
{
	int d = qn_mden(q);
	for (int k = 0; k < q->M; k++) { fb[k] = q->mi[k] % d; vfull[k] = q->mi[k] / d; }
}

static fam fam_new(int M)
{
	fam F = { M, 0, 1024, 4096, NULL, NULL };
	F.v = malloc((size_t)F.cap*M*sizeof(int)); F.key = malloc(F.hcap*sizeof(int));
	for (int z = 0; z < F.hcap; z++) F.key[z] = -1;
	return F;
}

/* U[r] for r = 1..R; returns the family size, or -1 above cap */
static long gcomom_unknowns(qnmodel* q, double* U, long cap)
{
	int M = q->M, R = q->R;
	fam F = fam_new(M);
	int* t = malloc(M*sizeof(int)); int* fb = malloc(M*sizeof(int)); int* vf = malloc(M*sizeof(int));
	split_mi(q, fb, vf);
	fam_add(&F, vf);
	for (int r = 1; r <= R; r++) U[r] = 0;
	long n = 0;
	for (int id = 0; id < F.n; id++) {
		if (F.n > cap) { n = -1; break; }
		int* v = F.v + (size_t)id*M, nt = 0, nrem = 0;
		for (int k = 0; k < M; k++) { if (v[k] > 0 || fb[k] > 0) nt++; if (v[k] > 0) nrem++; }
		int D = nrem >= R-1 ? 1 : (nt > 0 ? nt : 1);
		for (int r = 1; r <= R; r++) U[r] += nt * nckd(D + r - 1, r - 1);
		int nb = MIN(R-1, nrem), seen = 0;
		for (int k = M-1; k >= 0 && seen < nb; k--) {
			if (v[k] <= 0) continue;
			seen++;
			memcpy(t, F.v + (size_t)id*M, M*sizeof(int)); t[k]--;
			int any = 0; for (int z = 0; z < M; z++) if (t[z] || fb[z]) { any = 1; break; }
			if (any) fam_add(&F, t);
			v = F.v + (size_t)id*M;   /* fam_add may move F.v */
		}
	}
	if (n == 0) n = F.n;
	free(F.v); free(F.key); free(t); free(fb); free(vf);
	return n;
}

/* ---- gmom model tree (gmom/main.c gtree_new) on copy vectors: support d
 * = stations with a copy or a fractional base; a copy is removable unless
 * it is the last one of a one-station model without base; the GCEs remove
 * one copy of each of the last ntag = min(Bmax, #removable) of them; basis
 * level max(1, r - Bc), Bc = 1 for d = 1 else ntag.  Unit multiplicities
 * give the distinct-queue tree.  Returns sum_{r>=2} sum_models colA^3 +
 * N_r colA^2, INFINITY above cap models ---- */
static double gmom_work(qnmodel* q, int bmax, long cap)
{
	int M = q->M, R = q->R;
	if (R < 2) return INFINITY;
	fam F = fam_new(M);
	int* t = malloc(M*sizeof(int)); int* fb = malloc(M*sizeof(int)); int* vf = malloc(M*sizeof(int));
	int* rem = malloc(M*sizeof(int));
	split_mi(q, fb, vf);
	fam_add(&F, vf);
	double W = 0;
	for (int id = 0; id < F.n; id++) {
		if (F.n > cap) { W = INFINITY; break; }
		int* v = F.v + (size_t)id*M, d = 0, nrem = 0;
		for (int k = 0; k < M; k++) if (v[k] > 0 || fb[k] > 0) d++;
		for (int k = 0; k < M; k++)
			if (v[k] >= 1 && (d >= 2 || v[k] >= 2 || fb[k])) rem[nrem++] = k;
		int ntag = nrem < bmax ? nrem : bmax, Bc = (d == 1) ? 1 : ntag;
		for (int r = 2; r <= R; r++) {
			int l = r - Bc; if (l < 1) l = 1;
			double colA = r * nckd(d + l - 1, l);
			W += colA*colA*colA + (double)q->N[r-1] * colA*colA;
		}
		for (int c = 0; c < ntag; c++) {
			memcpy(t, F.v + (size_t)id*M, M*sizeof(int)); t[rem[nrem - ntag + c]]--;
			fam_add(&F, t);
		}
	}
	free(F.v); free(F.key); free(t); free(fb); free(vf); free(rem);
	return W;
}

/* ---- model features ---- */
typedef struct {
	int M, R, Ntot; double logb;
	bool ld, open, mi1, frac, zmixed, distinct_rows;
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
	f->frac = qn_mden(q) > 1;
	for (int k = 0; k < M; k++) { if (q->mi[k] != 1) f->mi1 = false; smi += q->mi[k]; }
	int nz = 0, zz = 0; for (int r = 0; r < R; r++) { if (mpz_sgn(q->Z[r])) nz++; else zz++; }
	f->zmixed = nz > 0 && zz > 0;
	if (f->ld || f->open) { f->logb = log(64.0); for (int s = 0; s < NSOLV; s++) f->W[s] = INFINITY; return; }
	/* bits of G (numerator, given the denominator bound) and the modulus
	 * that mommod/comommod draw primes for (util/gbound.h) */
	mpz_t* D = malloc((1 + q->R + q->M * q->R) * sizeof(mpz_t));
	long modbits = gb_modbits(q, D);
	double lg = gb_log2_G(q, q->N, -1);
	double bits = (isfinite(lg) && lg > 0 ? lg : 0) + (double)mpz_sizeinbase(D[0], 2);
	for (int i = 0; i < (1 + q->R + q->M * q->R); i++) mpz_clear(D[i]);
	free(D);
	f->logb = log(bits + 64.0);

	double ca = R * (double)smi; for (int r = 0; r < R; r++) ca *= (q->N[r] + 1.0);
	double mom = 0, binmom = 0, comom = 0, bincomom = 0;
	double gmom = gmom_work(q, 1, 400000), gmomb = R >= 3 ? gmom_work(q, R-1, 400000) : INFINITY;
	for (int r = 1; r <= R; r++) {
		int Nr = q->N[r-1];
		double dm = r * (nckd(M+r-1, r) + nckd(M+r-2, r-1));
		double dc = (M + 1) * nckd(M+r-1, r-1);
		double sm = dm*dm*M, ss = (double)Nr * dm * (M+1);
		mom += sm + ss;
		binmom += sm + ((Nr >= 2*dm*dm) ? dm*dm*dm*log2(Nr+1.0) : ss);
		comom += dc*dc*dc + (double)Nr*dc*dc;
		bincomom += dc*dc*dc + ((Nr - M - 1 >= 2*dc*dc) ? dc*dc*dc*log2(Nr+1.0) : (double)Nr*dc*dc);
	}
	double U[64]; double gc = 0;
	f->famsize = (R < 63) ? gcomom_unknowns(q, U, 400000) : -1;
	if (f->famsize > 0) for (int r = 1; r <= R; r++) gc += U[r] * (q->N[r-1] + (double)R*R) * R;
	double np = (modbits > 0 ? modbits : 2.0 * bits) / 60.0 + 1;   /* primes of 61 bits */
	f->W[S_CA] = ca; f->W[S_MOM] = mom; f->W[S_COMOM] = comom;
	f->W[S_GCOMOM] = f->famsize > 0 ? gc : INFINITY;
	f->W[S_GMOM] = gmom > 0 ? gmom : INFINITY;
	f->W[S_GMOMB] = gmomb > 0 ? gmomb : INFINITY;
	f->W[S_BINMOM] = binmom; f->W[S_BINCOMOM] = bincomom;
	f->W[S_MOMMOD] = mom * np; f->W[S_COMOMMOD] = comom * np;
}

static bool allowed(const feat* f, int s, bool dbg)
{
	if (!isfinite(f->W[s])) return false;
	/* fractional multiplicities: only the solvers that accept them */
	if (f->frac && s != S_MOM && s != S_COMOM && s != S_GCOMOM && s != S_GMOM && s != S_GMOMB) return false;
	switch (s) {
	case S_GMOM:   return f->R >= 2 && !dbg;
	case S_GMOMB:  return f->R >= 3 && !dbg;
	case S_BINMOM: return !f->zmixed;
	case S_MOMMOD: return !dbg;          /* no -d */
	default:       return true;
	}
}

static double predict(const feat* f, int s, int threads)
{
	/* fractional models: their own constants (larger coefficients, five
	 * candidates); a row left at zero by the fit falls back to calib.h */
	const calib_t* c = &CALIB[s];
	if (f->frac && (CALIB_FRAC[s].c0 != 0 || CALIB_FRAC[s].c1 != 0)) c = &CALIB_FRAC[s];
	double t = exp(c->c0 + c->c1 * log(fmax(f->W[s], 1.0)) + c->c2 * f->logb);
	/* primes are independent; the sub-models of one level in gcomom/gmomb
	 * scaled as about threads^0.6 on 8 threads (gmom's Bmax = 1 tree is a
	 * chain and does not scale) */
	if (threads > 1 && (s == S_MOMMOD || s == S_COMOMMOD)) t /= threads;
	if (threads > 1 && (s == S_GCOMOM || s == S_GMOMB)) t /= pow(threads, 0.6);
	return t + 0.003;
}

/* ---- running a solver ---- */
static char BINDIR[PATH_MAX];

static int GMOMB_B = 2;   /* R-1, set in main */

/* start candidate s; its stdout goes to outfd (-1: inherited), its stderr
 * to /dev/null if quiet */
static pid_t spawn(int s, char** flags, int nflags, const char* jarg, const char* model, int outfd, bool quiet)
{
	char path[PATH_MAX + 64], bb[16];
	snprintf(path, sizeof path, "%s/%s", BINDIR, BIN[s]);
	char* av[40]; int n = 0;
	av[n++] = (char*) BIN[s];
	for (int i = 0; i < nflags; i++) av[n++] = flags[i];
	if (s == S_GMOM || s == S_GMOMB) av[n++] = "-X";   /* never an approximate perturbation */
	if (s == S_GMOMB) { snprintf(bb, sizeof bb, "%d", GMOMB_B); av[n++] = "--bmax"; av[n++] = bb; }
	if (jarg && HAS_J[s]) { av[n++] = "-j"; av[n++] = (char*) jarg; }
	av[n++] = (char*) model; av[n] = NULL;
	fflush(stdout); fflush(stderr);
	pid_t p = fork();
	if (p == 0) {
		if (outfd >= 0) dup2(outfd, 1);
		if (quiet) { int dn = open("/dev/null", O_WRONLY); if (dn >= 0) dup2(dn, 2); }
		execv(path, av); fprintf(stderr, "automom: cannot run %s\n", path); _exit(127);
	}
	return p;
}

static int wait_rc(pid_t p) { int st; if (waitpid(p, &st, 0) < 0) return 128; return WIFEXITED(st) ? WEXITSTATUS(st) : 128; }

static void copy_out(int fd)
{
	char buf[65536]; ssize_t k;
	lseek(fd, 0, SEEK_SET);
	while ((k = read(fd, buf, sizeof buf)) > 0) fwrite(buf, 1, k, stdout);
	fflush(stdout);
}

static int run(int s, char** flags, int nflags, const char* jarg, const char* model)
{
	return wait_rc(spawn(s, flags, nflags, jarg, model, -1, false));
}

static int tmpfile_fd(void)
{
	char tmpl[] = "/tmp/automom.XXXXXX";
	int fd = mkstemp(tmpl);
	if (fd >= 0) unlink(tmpl);
	return fd;
}

/* buffered run: the output of a failed attempt is discarded */
static int run_capture(int s, char** flags, int nflags, const char* jarg, const char* model, bool quiet)
{
	int fd = tmpfile_fd();
	if (fd < 0) return run(s, flags, nflags, jarg, model);
	int rc = wait_rc(spawn(s, flags, nflags, jarg, model, fd, quiet));
	if (rc == 0) copy_out(fd);
	close(fd);
	return rc;
}

/* run a and b at once, keep the first to succeed and kill the other;
 * returns 0, or nonzero if both failed.  Both run with stderr discarded. */
static int race(int a, int b, char** flags, int nflags, const char* jarg, const char* model, bool explain)
{
	int fa = tmpfile_fd(), fb = tmpfile_fd();
	if (fa < 0 || fb < 0) { if (fa >= 0) close(fa); if (fb >= 0) close(fb); return -1; }
	pid_t pa = spawn(a, flags, nflags, jarg, model, fa, true);
	pid_t pb = spawn(b, flags, nflags, jarg, model, fb, true);
	int rc = 1;
	for (int left = 2; left > 0; left--) {
		int st; pid_t w = waitpid(-1, &st, 0);
		if (w != pa && w != pb) { left++; continue; }
		int ok = WIFEXITED(st) && WEXITSTATUS(st) == 0;
		bool isa = (w == pa);
		if (isa) pa = 0; else pb = 0;
		if (ok) {
			pid_t other = isa ? pb : pa;
			if (other > 0) { kill(other, SIGKILL); waitpid(other, NULL, 0); }
			if (explain) fprintf(stderr, "automom: race won by %s\n", NAME[isa ? a : b]);
			copy_out(isa ? fa : fb);
			rc = 0; break;
		}
	}
	close(fa); close(fb);
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
	printf("  --solver NAME    run NAME (one of ca mom comom gcomom gmom gmomb binmom bincomom mommod comommod;\n");
	printf("                   gmomb is gmom --bmax R-1)\n");
	printf("  --race           run the two best candidates at once when the second is predicted\n");
	printf("                   within 3x of the first, keep the first to finish (stderr discarded)\n");
	printf("  --features       print the cost features (used by automom/calibrate.py)\n");
}

int main(int argc, char** argv)
{
	char* flags[16]; int nf = 0;
	bool explain = false, dry = false, featmode = false, dbg = false, nc_only = true, racing = false;
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
		else if (!strcmp(a, "--race")) racing = true;
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

	qn_allow_rational_mi = 1;   /* fractional mi: routed to mom/comom/gcomom/gmom */
	qnmodel* q = readmodel(file);
	if (!q) { fprintf(stderr, "automom: cannot read %s\n", file); return 1; }
	nckinit(MAXNCKTABLE - 1, MAXNCKTABLE - 1);
	feat f; features(q, &f);
	GMOMB_B = q->R - 1 > 0 ? q->R - 1 : 1;
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
			/* outside the candidate table: a plain fork/exec */
			char path[PATH_MAX + 64]; snprintf(path, sizeof path, "%s/%s", BINDIR, seq[i]);
			char* av[24]; int na = 0; av[na++] = (char*) seq[i];
			for (int z = 0; z < nf; z++) av[na++] = flags[z];
			av[na++] = file; av[na] = NULL;
			int fd = seq[i+1] ? tmpfile_fd() : -1;
			fflush(stdout); fflush(stderr);
			pid_t p = fork();
			if (p == 0) { if (fd >= 0) dup2(fd, 1); execv(path, av); fprintf(stderr, "automom: cannot run %s\n", path); _exit(127); }
			int rc = wait_rc(p);
			if (rc == 0 && fd >= 0) copy_out(fd);
			if (fd >= 0) close(fd);
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
	int i0 = 0;
	if (racing && !force && n >= 2 && pt[ord[1]] < 3.0 * pt[ord[0]]) {
		if (explain) fprintf(stderr, "automom: racing %s and %s\n", NAME[ord[0]], NAME[ord[1]]);
		int rc = race(ord[0], ord[1], flags, nf, jarg, file, explain);
		if (rc == 0) return 0;
		if (rc > 0) { i0 = 2; if (i0 < n) fprintf(stderr, "automom: %s and %s failed or declined, trying %s\n", NAME[ord[0]], NAME[ord[1]], NAME[ord[2]]); }
	}
	for (int i = i0; i < n; i++) {
		int s = ord[i];
		int rc = (i < n-1) ? run_capture(s, flags, nf, jarg, file, s == S_GMOM || s == S_GMOMB || s == S_BINMOM)
		                   : run(s, flags, nf, jarg, file);
		if (rc == 0) return 0;
		if (i < n-1) fprintf(stderr, "automom: %s failed or declined (exit %d), trying %s\n", NAME[s], rc, NAME[ord[i+1]]);
	}
	return 1;
}
