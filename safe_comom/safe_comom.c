#include <stdio.h>
#include "safe_tiers.h"

/* safe_comom - CoMoM that always returns an exact answer, using the
 * hybrid MVA/CoMoM idea of Casale, "CoMoM" (IEEE TSE 2009), Appendix A:
 * when the CoMoM coefficient matrix is singular for a degenerate model,
 * interleave convolution with CoMoM so the degenerate directions are
 * handled by an MVA/convolution-style recursion.
 *
 * Realised as a tiered orchestrator over the exact solvers already in the
 * tree (bin/comom, then comom under permuted class orders, then bin/ca);
 * the tier machinery is shared with safe_mom and lives in
 * util/safe_tiers.c.
 *
 * The appendix's efficiency refinement -- peeling only the MINIMAL set of
 * degenerate classes by convolution while keeping CoMoM for the rest,
 * rather than falling all the way back to convolution in Tier 3 -- is a
 * further optimisation over this exact-but-tiered version; see README.
 */

static const safe_cfg CFG = {
	.name         = "safe_comom",
	.solver       = "comom",
	.solver_flags = "",     /* comom has no auto-perturbation to suppress */
	.probe_flags  = "-e",   /* errors instantly on a singular system */
	.fallback     = "ca",
	.permuted_flags = "-t -q",   /* both are indexed by class */
	.usage =
	  "  CoMoM that always returns the exact answer, via the hybrid\n"
	  "  MVA/CoMoM of TSE'09 Appendix A: plain CoMoM, then CoMoM with a\n"
	  "  permuted class order, then exact convolution for genuine\n"
	  "  loading degeneracies.\n"
};

int main(int argc, char** argv)
{
	return safe_tiers_main(&CFG, argc, argv);
}
