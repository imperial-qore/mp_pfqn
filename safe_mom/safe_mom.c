#include <stdio.h>
#include "safe_tiers.h"

/* safe_mom - MoM that always returns an exact answer.
 *
 * bin/mom degrades silently on a singular coefficient matrix: it warns and
 * re-solves a model perturbed at digit 20, so -g/-l/-t/-q return an
 * APPROXIMATE answer (with -e it refuses outright).  safe_mom keeps the
 * answer exact by applying the tiered hybrid MVA/MoM construction of
 * Casale, "CoMoM" (IEEE TSE 2009), Appendix A, to the MoM basis:
 *
 *   Tier 1  bin/mom -X            (no-perturb: exact or fail)
 *   Tier 2  bin/mom -X over the O(R) recursion-class choices
 *   Tier 3  bin/ca                (exact convolution)
 *
 * Every tier is exact, so the printed result is always the exact G/X/Q.
 * The tier machinery is shared with safe_comom (util/safe_tiers.c).
 *
 * Note on Tier 2: the loading-only oracle in util/pfqn_sing.h is exact for
 * the b=1 basis of gmom, not for mom's basis, so the class orders are
 * probed by running mom itself rather than screened by the oracle.
 *
 * Unlike safe_gmom, Tier 3 keeps mom's full domain: bin/ca handles
 * multiserver stations by station expansion, exactly.
 */

static const safe_cfg CFG = {
	.name         = "safe_mom",
	.solver       = "mom",
	.solver_flags = "-X",      /* never accept mom's perturbed fallback */
	.probe_flags  = "-e -X",   /* cheap exact-only probe: fails at once if singular */
	.fallback     = "ca",
	.permuted_flags = "-t -q",   /* both are indexed by class */
	.usage =
	  "  MoM that always returns the exact answer: plain MoM, then MoM with\n"
	  "  a permuted class order, then exact convolution for genuine loading\n"
	  "  degeneracies.  Never returns a perturbed (approximate) result.\n"
};

int main(int argc, char** argv)
{
	return safe_tiers_main(&CFG, argc, argv);
}
