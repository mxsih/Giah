---
title: Cavity DOE Project Notes
course: IEE 572
group: 19
tags: [IEE572, DOE, CFD, OpenFOAM, lid-driven-cavity]
---

# 2D lid-driven cavity: solver-settings DOE (IEE 572, Group 19)

## 1. The project in brief

### What Ghia did
- **Paper:** Ghia, Ghia & Shin (1982), *J. Comput. Phys.* 48(3), 387–411.
- It is a **numerical benchmark, not a physical experiment.** They solved the 2D incompressible Navier–Stokes equations for a square cavity whose top wall (the lid) moves at a constant speed.
- **Method:** vorticity–stream-function form, solved with a coupled strongly implicit multigrid method. Convective terms use first-order upwind with a deferred second-order correction, so the scheme is second-order accurate at convergence.
- **Grids and range:** uniform 129×129 grid for Re up to 3200, and 257×257 for Re = 5000, 7500 and 10,000.
- **How they checked their work:** Re = 100 was compared with earlier published numerical solutions, and Re = 400 was run on both 129×129 and 257×257 to show that 129×129 is adequate for moderate Re. As far as the paper shows, no such grid check was done at Re = 1000.
- **Their convergence definition affected the results.** The paper states that the criterion for convergence on the fine grids influenced the final solutions.
- Their velocity profiles became a standard reference, and most cavity codes are validated against them.

### What we will do
- Run the cavity in **OpenFOAM** (steady, laminar) under different solver settings.
- For each run, measure how far the simulated centerline velocities are from Ghia's tables. This single number is the response $E$.
- Use a designed experiment to find which settings matter, whether they interact, and which combination is best.

### Why
- Solver settings are usually chosen by habit. Their individual effects, and especially their interactions, are rarely quantified together.
- Knowing which settings drive accuracy tells us where compute time is worth spending.

### How DOE helps
- **Efficiency:** vary all four factors together (16 combinations) instead of one at a time.
- **Interactions:** shows whether one setting's effect depends on another (for example, a fine mesh may only pay off with the 2nd-order scheme).
- **Ranking and significance:** which settings matter and which are at noise level.
- **Blocking:** handles Reynolds number in a structured way, and lets us test whether effects change with it.
- **Confirmation:** predict the best settings from the model, then verify with extra runs.

### Design at a glance
| Item | Choice |
|---|---|
| Factor A, mesh | 32×32 (−) vs 128×128 (+) |
| Factor B, convection scheme | 1st-order upwind (−) vs 2nd-order (+) |
| Factor C, pressure–velocity coupling | SIMPLE (−) vs SIMPLEC (+) |
| Factor D, convergence tolerance | 1e-3 (−) vs 1e-6 (+) |
| Block | Reynolds number: 100, 400, 1000 (may become a 5th factor, see Step 4) |
| Runs | 16 treatments × 3 Reynolds numbers = 48 |
| Response $E$ | RMS error vs Ghia, as % of lid speed |

$$E=\frac{100}{U_{lid}}\sqrt{\frac{1}{30}\left[\sum_{i=1}^{15}\left(u_i^{sim}-u_i^{ref}\right)^2+\sum_{j=1}^{15}\left(v_j^{sim}-v_j^{ref}\right)^2\right]}$$

Lower is better. $E$ will not reach zero, because the Ghia tables carry their own discretization error.

---

## 2. The Ghia data we use

### What the data is
- The cavity is a unit square with lid speed 1, and Re = 1/viscosity.
- Ghia tabulates two velocity profiles through the geometric center of the cavity:
  - **Vertical centerline (x = 0.5): u vs y.** The lid drags fluid to the right at the top (u = 1 at y = 1). The return flow below is negative, with a minimum at lower y.
  - **Horizontal centerline (y = 0.5): v vs x.** The primary vortex turns clockwise, so v is positive on the left (flow rising) and negative on the right (flow sinking).
- Each profile has 17 points at non-uniform positions, chosen to include the local maxima and minima.
- Positions are exact multiples of 1/128. They are the nodes of Ghia's 129×129 grid, which is the same resolution as our 128×128 cell mesh.

### What we use and why
- **Re = 100, 400, 1000 only.** Ghia tabulates up to 10,000, but those are outside our design.
- **The 15 interior points per profile (30 in total).** The two end points per profile are boundary conditions (u = 1 at the lid and 0 at the bottom, v = 0 at both side walls). Every simulation matches them exactly, so they add zero error and would dilute $E$.
- Values below were read from the page images (Tables I and II of the paper) and saved in `ghia_Re100_400_1000.csv`. Spot-check them against the paper.

### Vertical centerline (x = 0.5): u-velocity vs y
| Ghia pt | y | Re = 100 | Re = 400 | Re = 1000 |
|---|---|---|---|---|
| 126 | 0.9766 | 0.84123 | 0.75837 | 0.65928 |
| 125 | 0.9688 | 0.78871 | 0.68439 | 0.57492 |
| 124 | 0.9609 | 0.73722 | 0.61756 | 0.51117 |
| 123 | 0.9531 | 0.68717 | 0.55892 | 0.46604 |
| 110 | 0.8516 | 0.23151 | 0.29093 | 0.33304 |
| 95 | 0.7344 | 0.00332 | 0.16256 | 0.18719 |
| 80 | 0.6172 | −0.13641 | 0.02135 | 0.05702 |
| 65 | 0.5000 | −0.20581 | −0.11477 | −0.06080 |
| 59 | 0.4531 | −0.21090 | −0.17119 | −0.10648 |
| 37 | 0.2813 | −0.15662 | −0.32726 | −0.27805 |
| 23 | 0.1719 | −0.10150 | −0.24299 | −0.38289 |
| 14 | 0.1016 | −0.06434 | −0.14612 | −0.29730 |
| 10 | 0.0703 | −0.04775 | −0.10338 | −0.22220 |
| 9 | 0.0625 | −0.04192 | −0.09266 | −0.20196 |
| 8 | 0.0547 | −0.03717 | −0.08186 | −0.18109 |

### Horizontal centerline (y = 0.5): v-velocity vs x
| Ghia pt | x | Re = 100 | Re = 400 | Re = 1000 |
|---|---|---|---|---|
| 125 | 0.9688 | −0.05906 | −0.12146 | −0.21388 |
| 124 | 0.9609 | −0.07391 | −0.15663 | −0.27669 |
| 123 | 0.9531 | −0.08864 | −0.19254 | −0.33714 |
| 122 | 0.9453 | −0.10313 | −0.22847 | −0.39188 |
| 117 | 0.9063 | −0.16914 | −0.23827 | −0.51550 |
| 111 | 0.8594 | −0.22445 | −0.44993 | −0.42665 |
| 104 | 0.8047 | −0.24533 | −0.38598 | −0.31966 |
| 65 | 0.5000 | 0.05454 | 0.05186 | 0.02526 |
| 31 | 0.2344 | 0.17527 | 0.30174 | 0.32235 |
| 30 | 0.2266 | 0.17507 | 0.30203 | 0.33075 |
| 21 | 0.1563 | 0.16077 | 0.28124 | 0.37095 |
| 13 | 0.0938 | 0.12317 | 0.22965 | 0.32627 |
| 11 | 0.0781 | 0.10890 | 0.20920 | 0.30353 |
| 10 | 0.0703 | 0.10091 | 0.19713 | 0.29012 |
| 9 | 0.0625 | 0.09233 | 0.18360 | 0.27485 |

### Caveats about the reference
- The tables contain at least two apparent misprints, both in columns we do not use (Re = 3200 at y = 0.4531 prints −0.86636, and Re = 10,000 at y = 0.5000 prints 0.03111 with no minus sign). Do not use those columns unchecked.
- Reference uncertainty is likely larger at Re = 1000 than at Re = 400, since the paper checked grid dependence only at Re = 400.
- OpenFOAM stores values at cell centers, while Ghia's points sit at cell vertices. Always interpolate the same way for every run, because a 32×32 mesh has coarse spacing near the lid and interpolation error mixes with discretization error.

---

## 3. Step-by-step: what has to be done and where DOE comes in

### Step 1. Pilot run and check against Ghia
- Set up the cavity (unit square, lid speed 1, viscosity = 1/Re).
- Run one baseline case: 128×128, 2nd-order, Re = 100.
- Sample both centerlines, interpolate to Ghia's points, and compute $E$.
- **What the result decides:**
  - If it matches Ghia closely (small $E$, profiles overlay), the setup is trustworthy and we continue.
  - If it does not, stop and fix the setup before any experimental run. Check boundary conditions, viscosity, sampling positions, interpolation, and the lid corners. Nothing built on a wrong baseline is usable.
- *DOE concept:* this is Chapter 1 step 1–3 (recognize the problem, choose factors and levels, choose the response), checked on one real case.

### Step 2. Check that the levels and response behave
- Run the extreme corners (coarse and loose vs fine and tight) at Re = 100 and Re = 1000.
- **Look for two things:**
  - Does $E$ differ clearly between levels? If the levels give nearly identical $E$, the factor ranges are too narrow and need widening.
  - Do any combinations fail to converge? Probe the hardest corner (32×32, 2nd-order, Re = 1000).
- Fix a convergence rule: residuals below the tolerance within a fixed iteration cap, plus a stable monitored quantity. Relaxation factors stay fixed, because changing them would change the treatment.
- *DOE concept:* choosing factor ranges (Chapter 1 step 2). Textbook section 6.5 advises spreading levels widely when the size of the error is unknown.

### Step 3. Exploratory check: do factor effects change with Reynolds number?
- Run a half fraction (8 runs, $2^{4-1}$) at Re = 100 and at Re = 1000, and compare the effect estimates.
- These runs are exploratory only and are not part of the final data.
- *DOE concept:* fractional factorials (Chapter 8) and interactions (Chapter 5). The grader asked us to check this before assuming an additive block effect.

### Step 4. Fix the design decisions
- **Reynolds number: block or factor?** If the effects change with Reynolds number, it becomes a fifth factor. The 48-run table is the same either way. Only the model changes.
- **Error term.** CFD is deterministic, so identical runs give identical results and there is no pure error. State this plainly and use one of these:
  - Pool the three-way and four-way interactions as error (sparsity of effects).
  - Use a normal probability plot of effects (Daniel's method, textbook section 6.5), combining the apparently negligible effects as error.
  - Use design projection: if some factors are negligible, the design collapses into a replicated factorial in the remaining factors (hidden replication).
  - Optionally, perturb the initial field to get small real variation at loose tolerance.
  - Be clear in the report that this error reflects model inadequacy, not random noise, so p-values are approximate.
- **Non-converged runs.** Record them, treat them as missing cells, and handle them by estimating the missing value from the model or by a documented reduced analysis. Report convergence success as a secondary outcome.
- *DOE concepts:* blocking (Chapters 4 and 7), unreplicated $2^k$ designs (Chapter 6).

### Step 5. Generate the design
- In JMP Custom Designer, create the $2^4$ factorial for each Reynolds number, with run order randomized within each block.
- *DOE concept:* randomization and blocking, two of the three basic principles (Chapter 1). Randomization matters most if runs share state, such as restart files.

### Step 6. Run the experiment
- Execute all 48 runs from the design table, compute $E$ for each, and record a convergence flag.
- The dataset is one table of 48 rows: run, block, A, B, C, D, $E$, converged.

### Step 7. Analyze in JMP
1. Estimate the effects and fit the model (Chapter 6).
2. Use a normal probability plot to separate real effects from noise.
3. Run the ANOVA with the chosen error term, including the Block × Factor interactions (Chapters 3 and 6).
4. Check residuals for normality and constant variance, and consider a log transform of $E$ (Chapter 3).
5. If some factors are negligible, project the design to get hidden replication.
6. Read the main-effects and interaction plots.

### Step 8. Confirm
- Run one to three extra cases at the predicted best settings and compare predicted vs actual $E$.
- Because the solver is deterministic, a large gap points to model inadequacy, not random noise.

### Step 9. Conclusions
- Which factors matter most, which interact, and the best setting combination.
- The cost–accuracy trade-off, such as whether the tight tolerance or the fine mesh is worth the extra compute.
- Scope: 2D laminar cavity, OpenFOAM, Re = 100–1000. Results should not be generalized beyond that.

### Where each DOE concept appears
| Concept (chapter) | Where used |
|---|---|
| 7-step procedure, principles (Ch 1) | Whole project; report structure |
| Hypothesis tests, ANOVA, assumptions, transformations (Ch 2–3) | Step 7 |
| Blocking, RCBD (Ch 4) | Reynolds number as block, Steps 4–5 |
| Factorial designs, interactions (Ch 5) | Design and analysis |
| $2^k$ designs, unreplicated analysis, normal plot (Ch 6) | Steps 4 and 7 |
| Blocking and confounding in $2^k$ (Ch 7) | Only if runs must be split into smaller blocks |
| Fractional factorials (Ch 8) | Step 3 exploratory check |
| Response surface (Ch 11) | Optional, after the main analysis |

---

## Report 2: what is needed and what is expected

- **Due:** Nov 9, 11:59 pm MST. Submit as "Report 2_Group 19", not zipped, by the group representative only.
- **Handout:** Report 2 is the design you plan to use, with factors, their roles (factor of interest or blocking), and the design.
- **Grader's comments on Report 1**, which the report must address:
  1. How the error term for statistical tests will be estimated when identical runs give identical results.
  2. Whether factor effects change with Reynolds number, checked before assuming an additive block effect.
  3. How runs that fail to converge will be handled.
- **Suggested contents:** short recap of Report 1 (with any level changes from the pilot), factor roles, the $2^4$ design with the treatment table, randomization and the run table, the error-term strategy, the Reynolds-number treatment, the non-convergence rule, and the planned analysis.
- **Do before writing:** Steps 1 to 5 above.
- **Final report:** due Dec 4, 10 pages maximum plus an appendix for code, with the dataset in a separate file.

## References
- U. Ghia, K. N. Ghia, C. T. Shin, "High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method," *J. Comput. Phys.* 48(3), 387–411 (1982).
- D. C. Montgomery, *Design and Analysis of Experiments*, Ch 6 (section 6.5) and Ch 7.
