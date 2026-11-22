#include <stdio.h>
#include "safe_tiers.h"

/* safe_promom - marginal queue-length probabilities, exact wherever the
 * class order can make them so.
 *
 * bin/promom inherits mom's coefficient matrix, so it is singular on exactly
 * the models bin/mom -X rejects and declines them rather than perturbing
 * (bin/procomom perturbs at digit 20 and returns an approximation).  This
 * wrapper applies the class-permutation tier of Casale, "CoMoM" (IEEE TSE
 * 2009), Appendix A:
 *
 *   Tier 1  bin/promom                (exact or fail)
 *   Tier 2  bin/promom over the O(R) recursion-class choices
 *   Tier 3  bin/camarg                (exact convolution marginals)
 *
 * Both tiers are exact, so a printed result is always the exact marginal.
 * The tier machinery is shared with safe_mom, safe_comom and safe_gmom
 * (util/safe_tiers.c).
 *
 * Tier 3 is bin/camarg, the convolution analogue for distributions: it has no
 * coefficient matrix, so it is exact on the degenerate models too.  It costs
 * the full population lattice, exponential in R, which is why it is the last
 * tier and not the first.  (This module originally had no Tier 3, bin/ca
 * returning means rather than distributions; camarg was written to close
 * that.)
 *
 * Output.  promom's -P (per station, sumN+1 probabilities) and -q (per
 * station, one mean) carry no class index: the marginal at a station does not
 * depend on the order the classes are processed in.  So permuted_flags is
 * empty and a Tier-2 result is printed exactly as it comes back, unlike the
 * -t/-q of the other wrappers.
 */

static const safe_cfg CFG = {
	.name           = "safe_promom",
	.solver         = "promom",
	.solver_flags   = "",       /* promom has no auto-perturbation to suppress */
	.probe_flags    = "-q",     /* fails at the first singular class transition */
	.fallback       = "camarg", /* exact convolution marginals, always exact */
	.permuted_flags = "",       /* -P and -q are per station, not per class */
	.flags_usage    = "[-P|-q]",
	.usage =
	  "  Marginal queue-length probabilities that are always exact: plain\n"
	  "  promom, then promom with a permuted class order, then exact\n"
	  "  convolution marginals (bin/camarg) for genuine degeneracies.\n"
	  "  Never returns a perturbed result.\n"
};

int main(int argc, char** argv)
{
	return safe_tiers_main(&CFG, argc, argv);
}
