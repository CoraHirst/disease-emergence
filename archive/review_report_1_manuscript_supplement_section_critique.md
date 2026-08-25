# Review 1: Section-by-Section Critique of Manuscript and Supplement

## Scope and approach
I reviewed:
- `products/manuscript/manuscript/Hirst_2026_mpox_Rxiv.pdf` (main manuscript)
- `products/manuscript/manuscript/Hirst_2026_mpox_SI_Rxiv.pdf` (supplement)
- Supporting project context from active analysis/code files for interpretability.

This report focuses on scientific clarity, internal logic, assumptions, and communication quality section by section.

## Executive summary
The manuscript presents a clear and policy-relevant framing: waning orthopox cross-immunity can increase both spillovers and within-human transmission, and reducing transmission (`R_e`) is far more impactful than reducing spillovers (`S`) in most realistic regimes. The central insight is strong.

Main weaknesses are:
1. Multiple places where notation and wording are inconsistent (notably `R_0` vs `R_e`, and emergence endpoint definitions).
2. Several strong assumptions are acknowledged but not quantified with uncertainty in the main text.
3. A few date/parameter statements are internally inconsistent and should be harmonized.
4. The supplement has equations/notation that are difficult to parse and should be cleaned for reproducibility.

---

## Main manuscript detailed critique

## Abstract
Strengths:
- Clear statement of motivation (post-eradication niche opening), mechanism (waning immunity), and intervention relevance.
- Identifies the policy asymmetry between prevention and post-emergence control.

Concerns and improvements:
- The abstract asserts that transmissibility increases matter much more than spillover increases, but does not briefly mention the modeling regime where this holds (near-critical `R_e` and low-to-moderate `P_E`). Add one sentence on domain-of-validity.
- "global threat" is plausible but broad; consider specifying the modeled criterion (e.g., emergence to sustained transmission and further adaptation).

## 1. Introduction
Strengths:
- Strong historical framing (smallpox eradication -> cessation of vaccination -> declining cross-immunity).
- Good linkage from classic eradication thresholds to contemporary emergence risk.

Concerns and improvements:
- There is a conceptual jump from historical mpox transmission patterns to current emergence risk; this would be stronger with one concise quantitative anchor in the Introduction (e.g., a reference estimate for current `R_e` range by context).
- The introduction promises empirical grounding but does not pre-state key data limitations (ascertainment, heterogeneous contact structure, clade-specific differences). Add these earlier.

## Results and Discussion (combined section)

### Clarification of reproductive numbers and emergence stages
Strengths:
- Clear distinction between `R_0` and `R_e` and how immunity versus adaptation affect them.
- The staged emergence schematic is intuitive and useful.

Concerns and improvements:
- There is an important notation inconsistency later in the manuscript/captions regarding whether emergence is defined by exceeding `R_e > 1`, or by reaching a high-transmissibility endpoint (`R_0` or `R_e` near 3). Define one primary endpoint early and maintain it uniformly.

### Waning immunity and `R_e` of mpox
Strengths:
- Transparent decomposition of immunity decline drivers.
- Useful simple formulation `R_e(t)=R_0(1-cF(t))`.

Concerns and improvements:
- The Figure 2 caption appears to mix vaccination stop years (mentions 1975 in one place but assumption tied to born-before-1980 elsewhere). This should be made internally consistent.
- The narrative states `R_0` is "in the range of 1 and likely slightly greater than one" from values `R_e≈0.3`, `c≈0.85`, `F≈0.9`; algebraically this implies `R_0≈1.28`. That may still be "slightly above 1" depending phrasing, but the text should show the implied value explicitly to avoid ambiguity.
- The same `c` is used to scale both transmission and spillover success. That assumption is convenient but biologically nontrivial; susceptibility, infectiousness, and exposure risk are not necessarily symmetric. This deserves explicit caveating in the main text (not only implicitly).

### Relationship between `R_e` and emergence (`P_E`)
Strengths:
- Core finding (steep dependence of `P_E` on `R_e` near 1) is compelling and biologically plausible.
- Discussion of mutation rate (`mu`) and fitness increment (`sigma`) effects is appropriate.

Concerns and improvements:
- The endpoint "emergence" appears to shift between "first supercritical" and "adaptation to high transmissibility". Keep separate symbols or terms (e.g., `P_{E,>1}` and `P_{E,high}`).
- Figure caption text appears to acknowledge a notation mismatch ("really this should be `R_e >= 3`"). This should be resolved in final prose, not left as a parenthetical correction.

### Pandemic potential (`P_P`) and heatmaps
Strengths:
- Good use of `P_P = 1-(1-P_E)^S` and low-probability approximation.
- Useful interpretation of policy-relevant contour bands.

Concerns and improvements:
- The statement that if expected emergence time is <10 years "we might expect the pathogen to already have emerged" is too strong without accounting for nonstationary parameters, detection lags, and surveillance differences.
- Consider adding uncertainty/credible bands around trajectories; current deterministic curves can overstate precision.

### Interventions and costs of delay
Strengths:
- Clear and policy-relevant comparison: prevention threshold vs post-emergence control burden.
- Correctly emphasizes broad utility of transmission-reducing interventions.

Concerns and improvements:
- The intervention section could more explicitly separate biological effect size from implementation feasibility/cost. The manuscript acknowledges this briefly; expanding that would improve policy realism.
- The vaccine discussion includes speculative pathways (e.g., pathology-only vaccines possibly increasing emergence risk). Keep but mark clearly as hypothesis unless supported by modeled sensitivity analysis.

## Discussion
Strengths:
- Thoughtful treatment of uncertainty (adaptive landscape, heterogeneity, behavior change, delay from mutation to recognition).
- Appropriately broadens implications to post-eradication emergence generally.

Concerns and improvements:
- One sentence appears semantically incorrect as written: referring to "mpox" as the only human disease eradicated; this should refer to smallpox.
- Some epidemiologic statements (case counts, clade-specific transmission characterization) should be date-stamped explicitly in text to avoid rapid obsolescence.
- The concluding argument would be stronger with a concise "what data would most reduce uncertainty" list (e.g., clade-specific transmission fitness landscape, cross-protection against transmission, spillover surveillance intensity).

## References and presentation quality
Strengths:
- Broad reference coverage spanning theory, immunology, and recent mpox literature.

Concerns and improvements:
- Several references are dynamic resources (e.g., CDC/WHO web pages) without access dates/versioning in the extracted text. Add retrieval dates.
- Minor typography/encoding and spacing artifacts appear in the PDF text extraction; verify final production PDF for readability and symbol rendering.

---

## Supplement detailed critique

## SI Section 1: Multi-type branching process formulation
Strengths:
- Appropriate framework and assumptions are mostly explicit.
- Good alignment with prior Antia et al.-style models.

Concerns and improvements:
- Some equations/notation are difficult to parse in the PDF text (especially the expression for `m-1` and the `R(m')` threshold notation). This section needs notation cleanup for unambiguous reproducibility.
- The numerical solution details are minimal: specify solver, initialization strategy, convergence criteria, and how alternate roots are handled.
- "Final type has arbitrarily high `R`" is acceptable computationally, but include a sensitivity statement showing results are insensitive beyond a threshold final `R`.

## SI Section 2: Robustness to `mu` and `sigma`
Strengths:
- Useful stress test direction; directly relevant to model claims.

Concerns and improvements:
- Robustness is shown visually but not summarized quantitatively. Add at least one metric (e.g., slope of `log P_E` vs `R_e` across parameter grid, or relative contribution decomposition).
- Figure annotation should more clearly define fixed parameters and which are varied in each panel.

## SI Section 3: Pre-adaptation pathway
Strengths:
- Important extension; acknowledges reservoir polymorphism/pre-adaptation possibility.

Concerns and improvements:
- The pre-adaptation mechanism is described qualitatively but not fully parameterized in text (frequency model in reservoir, fitness cost assumptions, and mapping to spillover `S`).
- Clarify when this pathway materially changes conclusions versus only shifting trajectory position.

## SI references
- Sufficient for the SI scope, but methods reproducibility would benefit from direct linkage to exact code functions/files in this repository.

---

## High-priority revision checklist
1. Harmonize emergence endpoint definitions and notation (`R_0` vs `R_e`, "supercritical" vs "high-transmissibility").
2. Resolve internal date inconsistency around vaccination cessation assumptions (1975 vs 1980 framing).
3. Explicitly show/calibrate implied `R_0` values from chosen `R_e`, `F`, and `c` examples.
4. Add uncertainty treatment (or explicit bounds) for key inferred quantities (`c`, `S_0`, trajectory timing).
5. Clean SI equation typography/notation and add computational solver details.
6. Correct wording errors that alter meaning (e.g., eradicated disease statement).
