#include <stdio.h>
#include "safe_tiers.h"

/* safe_gmom - generalized MoM that always returns an exact answer.
 *
 * bin/gmom degrades silently on a model that is singular under every class
 * order: it warns and re-solves a model perturbed at digit 20, so -g/-l/-t/-q
 * return an APPROXIMATE answer (with -e it refuses outright).  safe_gmom
 * keeps the answer exact by applying the tiered hybrid MVA/MoM construction
 * of Casale, "CoMoM" (IEEE TSE 2009), Appendix A, to gmom's b=1 basis:
 *
 *   Tier 1  bin/gmom -X          (no-perturb: exact or fail)
 *   Tier 2  bin/gmom -X over the O(R) recursion-class choices
 *   Tier 3  bin/ca               (exact convolution)
 *
 * Every tier is exact, so the printed result is always the exact G/X/Q.
 * The tier machinery is shared with safe_mom and safe_comom
 * (util/safe_tiers.c).
 *
 * Note on Tier 2: gmom already does this internally.  Its init-time oracle
 * (util/pfqn_sing.h) is exact for the b=1 basis, so gmom itself moves a
 * degenerate class to the recursion position and un-permutes X/Q on the way
 * out.  Tier 2 here is therefore expected never to fire; it is kept for
 * uniformity with safe_mom/safe_comom and as a backstop if the two class
 * searches ever disagree.
 *
 * Tier 3 also covers gmom's DOMAIN restrictions, not just its singularities:
 * gmom is limited to distinct single-server closed queues with R>=2, while
 * bin/ca handles multiserver stations exactly by station expansion and
 * single-class models directly.  A model outside gmom's domain fails Tier 1
 * and Tier 2 and is answered exactly by Tier 3.
 */

static const safe_cfg CFG = {
	.name         = "safe_gmom",
	.solver       = "gmom",
	.solver_flags = "-X",      /* never accept gmom's perturbed fallback */
	.probe_flags  = "-e -X",   /* cheap exact-only probe: fails at once if singular */
	.fallback     = "ca",
	.usage =
	  "  Generalized MoM (b=1) that always returns the exact answer: plain\n"
	  "  gmom, then gmom with a permuted class order, then exact convolution\n"
	  "  for genuine loading degeneracies and for models outside gmom's\n"
	  "  single-server domain.  Never returns a perturbed (approximate)\n"
	  "  result; use bin/gmom -p directly for that.\n"
};

int main(int argc, char** argv)
{
	return safe_tiers_main(&CFG, argc, argv);
}
