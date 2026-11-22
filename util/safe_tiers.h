#ifndef SAFE_TIERS_H
#define SAFE_TIERS_H

/* safe_tiers - shared tiered orchestrator for the "always exact" wrappers
 * (safe_comom, safe_mom), realising the hybrid MVA/MoM idea of Casale,
 * "CoMoM" (IEEE TSE 2009), Appendix A.
 *
 *   Tier 1  the primary exact solver on the model as given.
 *   Tier 2  the primary solver after permuting the class order.  A
 *           singularity that depends only on which class is the recursion
 *           class (the last one) is removed by processing a non-degenerate
 *           class last; O(R) orders are probed, not O(R!).
 *   Tier 3  a fallback exact solver (convolution) for genuine loading
 *           degeneracies that no class order removes.
 *
 * Every tier is exact, so the printed result is always the exact answer.
 *
 * G(N) is invariant to the class order, so -e/-g/-l need no un-permutation;
 * -t (per class) and -q (per station,class) from a permuted Tier-2 run are
 * mapped back to the original class order before printing.  Which flags need
 * that is per wrapper (`permuted_flags`): promom's -P and -q are per station
 * and carry no class index, so they are printed as they come back.
 */

typedef struct {
	const char* name;          /* wrapper name, used in diagnostics */
	const char* solver;        /* Tier 1/2 solver in the same bin/ dir */
	const char* solver_flags;  /* appended to every Tier 1/2 run ("" if none);
	                              use it to suppress the solver's own
	                              approximate fallback, e.g. "-X" for mom */
	const char* probe_flags;   /* cheap Tier-2 singularity probe, e.g. "-e" */
	const char* fallback;      /* Tier 3 exact solver, e.g. "ca"; NULL when
	                              no exact fallback exists for this output
	                              (marginal distributions have none) */
	const char* permuted_flags;/* space-separated output flags whose result is
	                              indexed BY CLASS and so must be mapped back
	                              after a Tier-2 run, e.g. "-t -q".  Flags not
	                              listed are class-order invariant and are
	                              printed verbatim.  NULL means "-t -q". */
	const char* flags_usage;   /* flag list for the usage line (may be NULL) */
	const char* usage;         /* extra usage text (may be NULL) */
} safe_cfg;

/* Full main() for a wrapper: parses [-e|-g|-l|-t|-q] model.qn, runs the
 * three tiers, prints the first exact result.  Returns a process exit
 * status. */
int safe_tiers_main(const safe_cfg* cfg, int argc, char** argv);

#endif
