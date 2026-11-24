/* Inner exact solver of oimom: bin/comom itself, its main.c compiled into
 * this translation unit with main renamed (as bincomom does), so that its
 * capture entry point comom_sg_direct can be called in-process.
 *
 * oimom_inner_solve returns the exported constants of the load-independent
 * network q in the util/singular.h layout,
 *   out[0] = G(N), out[1+r] = G(N-1_r), out[SG_IDX_GK(R,k,r)] = G^{+k}(N-1_r),
 * after comom's exact reduction (identical stations, proportional classes)
 * and, on a singular system, its tiered exact fallback.  0 on success. */
#define main comom_cli_main
#include "../comom/main.c"
#undef main

int OIMOM_THREADS = 1;   /* -j: primes of the interpolation tier solved concurrently */

int oimom_inner_solve(qnmodel* q, mpq_t* out)
{
	static int configured = 0;
	char tier[256];
	if (!configured) { comommod_config(OIMOM_THREADS, 61, 0, 0, 0); configured = 1; }
	sg_ops ops = { comom_sg_direct, comommod_interp, NULL, " per prime" };
	sg_qn_fix_nt(q);
	int rc = sg_solve_reduced(q, out, tier, sizeof tier, &ops, 1);
	if (rc == -1) rc = sg_solve_tiered(q, out, tier, sizeof tier, &ops, 0);
	if (getenv("OIMOM_DEBUG")) {
		int r; fprintf(stderr, "inner R=%d M=%d N=(", q->R, q->M);
		for (r = 0; r < q->R; r++) fprintf(stderr, "%d%s", q->N[r], r < q->R-1 ? "," : "");
		fprintf(stderr, ") route: %s\n", tier);
	}
	return rc;
}
