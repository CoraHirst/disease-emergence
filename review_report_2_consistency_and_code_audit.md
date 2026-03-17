# Review 2: Consistency Check and Code Audit (Manuscript, Supplement, Code Base)

## Scope
I performed a second-pass audit focused on:
- Internal consistency across manuscript, supplement, and code.
- Potential implementation errors in the active analysis pipeline.
- Evidence of runtime issues from rendered outputs already committed in the repo.

Reviewed files include:
- Manuscript PDFs in `products/manuscript/manuscript/`
- Core code in `code/functions/*.R`
- Analysis notebooks in `code/data analysis/congo/*.qmd`
- Figure-generation notebook in `code/project/paper figures/paper_figures.qmd`
- Rendered evidence in `code/project/paper figures/paper_figures.html`
- Project READMEs.

## Limitations
- I could not re-run R code in this environment (R executable unavailable), so this is a static + rendered-output audit.
- PDF text extraction introduces minor encoding noise; findings below prioritize issues that remain clear despite extraction artifacts.

---

## Findings (ordered by severity)

## Critical

### C1. `waning_population_immunity_1.qmd` has an undefined object in executable code.
- Evidence: [`code/data analysis/congo/waning_population_immunity_1.qmd:166`](code/data analysis/congo/waning_population_immunity_1.qmd:166) uses `m` in `age_structure = data.frame(... m[, ...])`.
- No prior definition of `m` appears in this file.
- Impact: the age-structure chunk cannot run in a clean session, breaking notebook reproducibility.
- Recommendation: either compute `m` in this notebook or remove/guard the chunk.

## High

### H1. Contradictory probability formula in figure notebook prose.
- Evidence: [`code/project/paper figures/paper_figures.qmd:180`](code/project/paper figures/paper_figures.qmd:180) states `P = p(1-p)^{n-1}`.
- But code function uses `P = 1-(1-p)^n` in [`code/functions/functions.R:127`](code/functions/functions.R:127), and manuscript Eq. (4) also uses that form.
- Impact: conceptual inconsistency between narrative and computation; can mislead readers and reviewers.
- Recommendation: replace text at line 180 with the correct cumulative probability expression.

### H2. Internal date inconsistency for vaccination cessation assumptions.
- Evidence in manuscript text: Figure 2 caption mentions vaccination stopped in 1975, while nearby text and assumptions refer to cessation around 1980 and immunity as "born prior to 1980."
- Code assumptions are tied to 1980 (e.g., [`code/data analysis/congo/matrix_model_projections_2.qmd:174`](code/data analysis/congo/matrix_model_projections_2.qmd:174)).
- Impact: weakens parameter interpretation (`F(t)`, `c` calibration) and creates avoidable reviewer confusion.
- Recommendation: pick one operational definition (e.g., 1980) and consistently apply/report it.

### H3. Eradication wording error changes meaning.
- Evidence in manuscript discussion text: statement reads as if mpox is the eradicated disease (should be smallpox).
- Impact: factual inaccuracy in a key interpretive section.
- Recommendation: correct wording explicitly.

### H4. Emergence endpoint notation is inconsistent (`R_0` vs `R_e`, threshold vs high-transmissibility endpoint).
- Evidence:
  - Manuscript/captions mix "evolve past 1" with "reach 3" and includes note-like text indicating notation uncertainty.
  - Code computes mutation steps using `supercritical_R = 3` then sets final type to very high `R_adapted` (e.g., [`code/functions/functions.R:76-90`](code/functions/functions.R:76)).
- Impact: ambiguity in what exactly probabilities represent (`first supercritical` vs `high-R persistence`).
- Recommendation: define and label separate outcomes explicitly throughout.

## Medium

### M1. Inconsistent nonlinear solver initialization across closely related functions.
- Evidence:
  - [`code/functions/functions.R:31`](code/functions/functions.R:31) sets `xstart[1] = 1` in `pEmergence_deltaR`.
  - [`code/functions/functions.R:105`](code/functions/functions.R:105) sets `xstart[1] = 0` in `pEmergence_supercrit_deltaR`.
- Impact: root-finding can converge to different fixed points depending on initialization; risk of subtle probability bias.
- Recommendation: standardize initialization and verify against analytic/simulation checks for representative parameter points.

### M2. Rendered figure notebook shows repeated non-finite warnings (`NaN`, dropped rows).
- Evidence in rendered output:
  - [`code/project/paper figures/paper_figures.html:14125`](code/project/paper figures/paper_figures.html:14125) and many later lines show `NaNs produced`.
  - Multiple `Removed ... rows containing non-finite` warnings (e.g., [`...:39300`](code/project/paper figures/paper_figures.html:39300)).
- Likely source: `log10(P)` when `P=0` (or numerically underflowed), plus intentional `NA` endpoints in `geom_segment` arrows.
- Impact: warnings are expected in parts but currently noisy and unquantified; can hide real issues.
- Recommendation: clamp probabilities (e.g., `pmax(P, .Machine$double.xmin)` before log), and avoid constructing NA endpoints inside mapped aesthetics.

### M3. Off-by-one visual split between historical and projected periods.
- Evidence:
  - Historical styling uses `Year <= 2024` in figure code (e.g., [`paper_figures.qmd:67`](code/project/paper figures/paper_figures.qmd:67), [`:125`](code/project/paper figures/paper_figures.qmd:125)).
  - Projection styling uses `Year > 2023`.
- Impact: year 2024 appears in both conceptual groups depending on context; mismatch with manuscript statement "after 2023".
- Recommendation: use one consistent boundary (typically `<=2023` historical, `>=2024` projection).

### M4. Parameter sweep includes biologically/ numerically awkward boundary (`R_wt = 0`).
- Evidence: [`paper_figures.qmd:193`](code/project/paper figures/paper_figures.qmd:193) sets `R_wts = seq(0,1,0.01)` and function internally uses `log10(R_wt + 1e-8)`.
- Impact: extreme mutation-step counts and unstable tails near zero; contributes to warning-heavy output.
- Recommendation: start sweeps at a positive lower bound (e.g., 0.05 or 0.1) unless zero is explicitly justified.

### M5. Verbose solver tracing is enabled globally in core functions.
- Evidence: repeated `control=list(trace=1, stepmax=2)` in core routines (e.g., [`functions.R:35`](code/functions/functions.R:35), [`:57`](code/functions/functions.R:57), [`:71`](code/functions/functions.R:71), [`:109`](code/functions/functions.R:109)).
- Impact: massive stdout noise in rendered artifacts and slower runs.
- Recommendation: default `trace=0` (or parameterize verbosity).

## Low

### L1. Documentation path and spelling inconsistencies reduce reproducibility confidence.
- Evidence:
  - [`code/data analysis/README.md:14`](code/data analysis/README.md:14) path typo (`congo_demographics/1980-2023.xlsx` vs actual filename with underscore).
  - [`code/data analysis/README.md:23-24`](code/data analysis/README.md:23) includes misspellings (`proprotion`, `suriving`, `poulation_matrix.csv`).
- Impact: onboarding friction and avoidable confusion.
- Recommendation: clean README paths/spelling.

### L2. Plot labeling typo in age-structure panels.
- Evidence: `p_2023` chart title still says `Age Structure, 2000` in both data-analysis notebooks (e.g., [`waning_population_immunity_1.qmd:190`](code/data analysis/congo/waning_population_immunity_1.qmd:190), [`matrix_model_projections_2.qmd:250`](code/data analysis/congo/matrix_model_projections_2.qmd:250)).
- Impact: minor presentation error.
- Recommendation: correct label to 2023.

---

## Cross-document consistency matrix

1. Manuscript vs code on spillover/emergence math:
- Status: Mostly consistent in actual equations and implementation (`1-(1-p)^n`).
- Exception: prose error in `paper_figures.qmd` line 180.

2. Manuscript vs code on time anchoring:
- Status: Partially inconsistent (1975 mention vs 1980-based implementation).

3. Manuscript/SI vs code on emergence endpoint:
- Status: Conceptually close but not cleanly harmonized in notation and wording.

4. SI formalism vs implementation details:
- Status: Implemented model structure aligns broadly, but SI lacks sufficient numerical method detail and has notation readability issues.

---

## Recommended remediation order
1. Fix C1 and H1 immediately (reproducibility + mathematical clarity).
2. Harmonize date assumptions and endpoint notation (H2/H4).
3. Resolve solver/NaN handling issues (M1/M2/M4/M5) for stable reruns and cleaner figures.
4. Clean documentation and labels (L1/L2).

---

## Bottom line
The core modeling direction is coherent and valuable, but the current package has several high-visibility inconsistencies and one hard execution bug in an analysis notebook. Addressing the listed items would materially improve reviewer confidence in both scientific and computational reproducibility.
