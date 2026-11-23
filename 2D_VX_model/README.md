# 2D Vertex Model of Epithelial Wound Closure

A MATLAB implementation of a two-dimensional vertex model for epithelial tissue mechanics and wound closure. The tissue is a periodic network of polygonal cells whose vertices move under passive cell-mechanical forces and active wound-closing forces. Topological events (neighbour exchanges and wound-margin rearrangements) are selected using the local mechanical state and accepted through an energy criterion. The wound boundary is tracked directly as an ordered polygon throughout the simulation.

The folder also contains an independent validation of the implemented mechanics against the Farhadifar et al. (2007) soft/hexagonal phase diagram.

---

## Contents

1. [Overview](#1-overview)
2. [Model formulation](#2-model-formulation)
3. [Final results](#3-final-results)
4. [Getting started](#4-getting-started)
5. [Reproducing the final results](#5-reproducing-the-final-results)
6. [Other simulation scripts](#6-other-simulation-scripts)
7. [Parameters](#7-parameters)
8. [Output files](#8-output-files)
9. [File reference](#9-file-reference)
10. [Main data structures](#10-main-data-structures)
11. [References](#11-references)

---

## 1. Overview

| Item | Value |
| :--- | :--- |
| Domain | Periodic square box, `Lx = Ly = 15` |
| Initial tissue | Periodic Voronoi tessellation of 225 cells |
| Wound | All cells with centres inside a circle of radius 2 at the box centre are removed, leaving 213 cells |
| Tissue regime | Fluid-like, preferred shape index `p0 = 4` |
| Passive mechanics | Area elasticity and perimeter elasticity of every cell |
| Active mechanics | Purse-string line tension on the wound edges and inward lamellipodial crawling of the margin cells |
| Recruitment | Purse string, crawling and margin stiffening are recruited smoothly in time on a decaying clock `Kw(t)` |
| Topological events | Bulk T1 transitions, wound-margin handoffs, wound-corner trims (cell division and wound sealing available as options) |
| Event rule | Bell-type force-dependent rate compared with a deterministic threshold, followed by an energy-acceptance test |
| Dynamics | Overdamped vertex motion, explicit Euler integration, `dt = 0.01`, 1000 steps |
| Wound measurement | Exact shoelace area of the ordered wound loop at every timestep |
| Validation | Force/energy-gradient consistency test and reproduction of the Farhadifar soft/hexagonal phase boundary |

The main simulation is `runForceClosureStrategy.m`. The phase-diagram validation is `phaseDiagramValidation.m`.

---

## 2. Model formulation

All quantities are in normalised model units (preferred cell area `A0 = 1`).

### 2.1 Tissue representation

Each cell is an ordered, anticlockwise list of vertex indices; each vertex stores a 2D position. Opposite sides of the box are identified (periodic boundary conditions), so every distance, area and perimeter is computed with periodic corrections (`modifyVerticesForPeriodicity.m`). The initial tissue is a periodic Voronoi tessellation generated from randomly placed seeds with a minimum spacing (`genDataClassical.m`). All driver scripts call `rng('default')`, so a given parameter set always produces the same tissue and the same trajectory.

### 2.2 Cell energy and passive forces

The mechanical energy of the cells is

```
E_cells = sum_c [ k_c (A_c - 1)^2 + (1/r_c) (P_c - p0)^2 ]
```

where `A_c` and `P_c` are the area and perimeter of cell `c`, `p0` is the preferred shape index (equal to the preferred perimeter since `A0 = 1`), `k_c` is the area coefficient and `1/r_c` is the perimeter coefficient. Bulk values are `k_a = 1`, `r_s = 0.5`, `p0 = 4`.

The force that cell `c` exerts on vertex `i` is the negative energy gradient,

```
F_i^(c) = -2 k_c (A_c - 1) dA_c/dr_i - (2/r_c) (P_c - p0) dP_c/dr_i
```

with analytical derivatives from `getAreaDerivative.m` and `getPerimeterDerivative.m`. Contributions of all cells sharing a vertex are summed in `getVertexForcesClassical.m`.

### 2.3 Active wound forces

- **Purse string.** Every free wound edge carries a line tension `lambda`, with energy `E_ps = lambda * sum(L_ij)` over the wound edges. Each wound edge pulls its two end vertices towards each other with force `lambda * e_ij`; along the curved boundary this contracts the wound.
- **Crawling.** Each wound-margin cell receives an inward force of magnitude `F_crawl`, directed from the cell centre towards the wound centre and shared equally among the cell's vertices.
- **Margin stiffening.** The area coefficient and perimeter contractility of the margin cells are multiplied by factors that rise from 1 towards `ka_wound_factor = 2` and `contractility_wound = 3`.

### 2.4 Recruitment of the active forces

A decaying clock `Kw(t) = Kw0 exp(-t/tau_Kw)` defines the recruited fraction `g = 1 - Kw/Kw0`. Each active amplitude relaxes towards its target at rate `k_recruit`:

```
lambda      <- lambda + dt k_recruit (lambda_purse_string * g - lambda)
ka factor   <- relaxes towards 1 + (ka_wound_factor - 1) g
contr. fac. <- relaxes towards 1 + (contractility_wound - 1) g
F_crawl     =  crawl_force0 * max(0, 1 - lambda / lambda_purse_string)
```

so crawling is strongest just after wounding and hands over to the purse string as it assembles. With `lambda_purse_string = 4`, `crawl_force0 = 3`, `k_recruit = 1` and `tau_Kw = 2`, the margin factors reach 1.987 and 2.973 and the tension reaches 3.947 by step 1000.

### 2.5 Optional wound-area resistance

When `param.enableWoundResistance = true`, `Kw(t)` also acts as an elastic stiffness of the wound opening, `E_w = Kw(t) (A_w - A_w0)^2`, with the matching force on the wound-loop vertices. In the default configuration the clock is used for recruitment only, so closure is driven entirely by the purse string and crawling.

### 2.6 Overdamped dynamics

Vertices move according to `eta dr_i/dt = F_i`, integrated by explicit Euler (`updateVertexPositions.m`):

```
r_i(n+1) = r_i(n) + (dt/eta) F_i(n),   dt = 0.01, eta = 1
```

A second integrator, `updateVertexPositionsMonotonic.m`, halves the step until the mechanical energy does not increase; it is selected with `param.enforceMonotonicEnergy = true` (that mode also fixes the active coefficients and switches crawling off).

### 2.7 Topological events

All events follow the same three-stage rule:

1. **Geometric eligibility** (a short edge, a crowded margin, a small gap).
2. **Force-dependent Bell rate.** With a local mechanical drive `F`,
   ```
   k_off = k_off0 exp(F / f_beta_eff),   P_event = 1 - exp(-dt k_off)
   ```
   and the event is a candidate when `P_event` exceeds a fixed threshold. The threshold comparison makes event selection deterministic and the runs exactly reproducible.
3. **Energy acceptance.** The event is built as a trial configuration and its energy change `dE` is evaluated with `getTissueEnergyClassical.m`. With the Metropolis temperature `T_eff = 0`, an event is accepted when `dE < 0`; a positive `T_eff` additionally accepts energy-increasing events with probability `exp(-dE/T_eff)`.

**Bulk T1 transitions** (`checkT1transitions.m`, `performT1swap.m`). An edge shorter than `T1_TOL_bulk = 0.05` (or `T1_TOL_wound = 0.15` for margin cells) is a candidate. The drive is the edge tension per unit length, with `f_beta_eff = max(f_beta_T1 - (TOL - L), 0.01)`. One T1 is accepted per step, and its vertices then rest for `maxT1relaxsteps = 30` steps.

**Wound-margin handoffs and corner trims** (`checkWoundIntercalations.m`, `performWoundMarginT1.m`, `trimWoundCorner.m`). The wound is treated as the fourth cell of a T1 tetrad. In a *handoff*, a margin cell transfers its wound contact to its two neighbours and withdraws into the bulk while remaining intact; the new edge starts at length `woundT1NewEdgeLength = 0.01` and grows under the forces. A *corner trim* removes a redundant corner vertex where two wound edges meet. The drive is the compressive area pressure of the margin cell. Handoffs are paced by the measured contraction of the wound perimeter,

```
N_target = max(5, ceil(N_b0 * min(1, P_w / P_w0)))
```

and a handoff is allowed while the margin holds more than `N_target` cells. The Bell threshold moves from `0.0125` towards `0.009` as the perimeter contracts. Handoffs are tried before trims. The energy used by the margin test includes `woundBoundaryCellEnergy = 5` per boundary cell, which carries no positional force and favours turnover of the margin.

**Cell division** (optional, `checkCellDivision.m`, `performCellDivision.m`). A cell larger than `Area_div_thresh` is a candidate, driven by its area pressure `2 k (A - 1)`; it divides through its centroid along its short axis (`getCellPrincipalAxis.m`).

**Wound sealing** (optional, `checkWoundSealing.m`, `performWoundSeal.m`). Two opposing margin vertices closer than `L_seal_thresh` are merged at their midpoint once the margin has at most `maxBoundaryCellsForSealing` cells.

### 2.8 Wound measurement and diagnostics

After every timestep the free edges are found (`findWoundEdgeCells.m`, hash map of edges), ordered into a closed loop (`getOrderedWoundLoop.m`) and the enclosed area is computed exactly by the shoelace formula

```
A_w = 1/2 sum_i (x_i y_(i+1) - y_i x_(i+1))
```

`getWoundDiagnostics.m` also records the box-area deficit `Lx Ly - sum(A_c)`, the wound perimeter, the boundary cell/edge/vertex counts, vertex degrees, the smallest cell area and edge length, and a connectivity check from `validateTissueTopology.m`.

### 2.9 Order of operations in one timestep

1. Update the recruitment clock and the active amplitudes.
2. Compute vertex forces and advance the vertices.
3. Update cell areas and perimeters; compute edge forces (`getEdgeForces.m`).
4. Wound-margin handoffs/trims (recompute forces if one is accepted).
5. Bulk T1 transitions (refresh the wound margin if one is accepted).
6. Cell division and wound sealing (when enabled).
7. Wound diagnostics; snapshot every 50 steps.

---

## 3. Final results

### 3.1 Wound-closure simulation

Results folder: `results/endsem_btp_results/`

| Quantity | Value |
| :--- | :--- |
| Steps / model time | 1000 / 10 |
| Cells | 213 throughout |
| Initial wound-loop area | 12.33 |
| Peak loop area | 106.6% of initial, at step 117 |
| 50% / 25% / 1% of initial area | steps 350 / 427 / 515 |
| Final wound-loop area | 2.20 x 10^-3 (0.018% of initial; 99.98% reduction) |
| Wound-margin cells | 15 -> 5 |
| Accepted bulk T1 transitions | 17 |
| Wound-margin handoffs / corner trims | 10 / 9 |
| Mechanical energy | 105.2 -> 27.6 |
| Final margin factors (area / contractility) | 1.987 / 2.973 |
| Final purse-string tension | 3.947 |

### 3.2 Phase-diagram validation (MATLAB)

Results folder: `results/phase_diagram_validation_2026-11-19_185724/`

Matching the cell energy to the Farhadifar form `K/2 (A - A0)^2 + Gamma/2 P^2 + Lambda * sum(L)` gives `K = 2 ka`, `Gamma = 2/rstiff`, `Lambda = -4 p0/rstiff`, so that

```
Gbar = 1/(ka rstiff),   Lbar = -2 p0 Gbar
```

and the soft/hexagonal boundary is the line `p0 = sqrt(8 sqrt(3)) = 3.7224`.

| Check | Result |
| :--- | :--- |
| Analytical force vs. numerical energy gradient | relative difference 1.837 x 10^-10 |
| Points relaxed | 55 (35 hexagonal, 20 soft), every point in its theoretical phase |
| Boundary along the line (three Gbar rows) | hexagonal at p0 = 3.720, soft at p0 = 3.725, bracketing 3.7224 |
| Soft networks | relaxed shape index settles at p0 |
| Hexagonal networks | residual energy per cell matches `Gbar (3.7224 - p0)^2` |

### 3.3 Independent NumPy cross-check

Results folder: `results/phase_diagram_validation_numpy_port_2026-11-19/`

A NumPy reimplementation of the bulk mechanics (`vxport.py`) relaxes 285 configurations and reproduces the same phase boundary; its forces match the numerical energy gradient to about 10^-9. All points shared with the MATLAB sweep have the same phase. See the `README.txt` in that folder.

---

## 4. Getting started

**Requirements:** MATLAB R2020a or newer (base MATLAB, no toolboxes). Python 3 with NumPy, pandas and Matplotlib is needed only for the optional NumPy cross-check.

**Working folder:** run every script from inside `2D_VX_model/`:

```matlab
cd('/path/to/woundclosure/2D_VX_model')
```

The drivers add `functions2D/` to the path themselves. Several helper names exist both in this folder and in `functions2D/`; MATLAB gives files in the current folder priority over the path, so running from `2D_VX_model/` uses the versions in this folder, which are the ones used for all results. `functions2D/createWound.m` and `functions2D/findWoundEdgeCells.m` exist only in `functions2D/`.

---

## 5. Reproducing the final results

### 5.1 Wound-closure simulation

```matlab
cd('/path/to/woundclosure/2D_VX_model')
setenv('WOUND_RUN_TAG', 'crawling_active_1000_steps');   % optional: only sets the folder name
runForceClosureStrategy
```

This writes `results/force_closure_gradual_intercalation_crawling_active_1000_steps_<timestamp>/` with the same numbers as Section 3.1. Without the tag, the folder is `results/force_closure_gradual_intercalation_<timestamp>/`.

Optional environment variables read by `runForceClosureStrategy.m`:

| Variable | Effect |
| :--- | :--- |
| `WOUND_NSTEPS` | number of timesteps (default 1000) |
| `WOUND_RUN_TAG` | text inserted in the output folder name |
| `WOUND_CONFIG_JSON` | path to a JSON file whose fields override `param` fields (unknown names stop the run) |
| `WOUND_SAVE_VIDEO` | `1` to record an MPEG-4 video with one frame per timestep |

Clear a variable with `setenv('WOUND_RUN_TAG', '')`.

### 5.2 Contractility and energy post-processing

`plot_contractility_annotations.m` (inside the final results folder) reads the `workspace.mat` next to it and writes the contractility and tension plots, the annotated energy plot, the per-step CSV and the notes file. To apply it to a new run, copy it into that run's folder and run it:

```matlab
run('results/<run folder>/plot_contractility_annotations.m')
```

### 5.3 Phase-diagram validation

```matlab
cd('/path/to/woundclosure/2D_VX_model')
phaseDiagramValidation
```

`runMode = 'quick'` (default, 55 points, about 25 minutes) reproduces Section 3.2; `runMode = 'full'` runs a denser sweep. Output: `results/phase_diagram_validation_<timestamp>/`.

### 5.4 NumPy cross-check (optional)

```bash
cd results/phase_diagram_validation_numpy_port_2026-11-19
python3 sweep.py        # writes sweep_results.csv and cfg_*.npz
python3 plot_phase.py   # writes phase_diagram_validation.png
```

---

## 6. Other simulation scripts

| Script | Purpose | Output folder |
| :--- | :--- | :--- |
| `runEnergyParameterTrial.m` | Function form of the main simulation: `runEnergyParameterTrial(configPath, tag)` applies the JSON overrides in `configPath`, runs, and also saves `metrics.json`, `parameters.json` and `energy_budget.csv` | `results/force_closure_gradual_intercalation_<tag>_<timestamp>/` |
| `continueEnergyTrial.m` | `continueEnergyTrial(workspacePath, extraSteps)` loads a saved run and continues it for `extraSteps` more steps | `results/force_closure_gradual_intercalation_extended_<timestamp>/` |
| `runMain.m` | Full closure run with video, snapshots and closure metrics (division and sealing follow the parameter file) | `testing/` |
| `runForceBasedTransitions.m` | Wound run instrumented for edge forces: force-vs-length scatter at T1 events, Trepat-style radial/azimuthal profiles, tension maps, force histograms | `results/wound_contractility/` |
| `runEdgeForceAnalysis.m` | Wound run followed by edge-force decomposition (perimeter tension, purse-string tension, area pressure) and force-vector plots | `results/edge_forces_new/` |
| `myMainClassical.m` | Wound run with a live plot of the vertex forces every `plotForces_Interval` steps | `results/myMainClassical_*` |
| `myMainClassicalWoundOnly.m` | Builds a wounded 30 x 30 tissue and plots its initial state | none |
| `mainClassical.m` | Classical vertex model without a wound (T1 transitions only) | `results/mainClassical_*` |

---

## 7. Parameters

All defaults are in `SingleCellContractionParameters.m`, which needs `param.Lx`, `param.Ly` and `Case` (`'fluid'`: `p0 = 4`; `'solid'`: `p0 = 2`) to be set first. Drivers then override a few fields (`runForceClosureStrategy.m`: `Nsteps = 1000`, `plotForces = false`, circular wound of radius 2 at the box centre).

**Cell mechanics and time stepping**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `ka` | 1.0 | bulk area coefficient |
| `rstiff` | 0.5 | bulk perimeter parameter (perimeter coefficient `1/rstiff`) |
| `p0` | 4.0 (fluid) / 2 (solid) | preferred shape index |
| `eta` | 1 | vertex drag coefficient |
| `deltat` | 0.01 | timestep |
| `Nsteps` | 3000 (1000 in the main driver) | number of timesteps |
| `isBoundaryFixed` | 0 | 0 = fully periodic box |

**Active forces and recruitment**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `lambda_purse_string` | 4.0 | target purse-string line tension |
| `t_ramp_purse_string` | 50 | linear ramp length, used when no recruitment state is present |
| `enableWoundCrawling` | true | crawling on/off |
| `crawl_force0` | 3.0 | initial crawling amplitude per margin cell |
| `tau_decay_crawl` | 100 | exponential crawl decay, used when no recruitment state is present |
| `ka_wound_factor` | 2.0 | target area-coefficient multiplier for margin cells |
| `contractility_wound` | 3.0 | target contractility multiplier for margin cells |
| `Kw0` | 2.0 | initial value of the clock `Kw` |
| `tau_Kw` | 2 | decay time of `Kw` |
| `k_recruit` | 1.0 | recruitment rate |
| `enableWoundResistance` | false | true: `Kw` also acts as wound-area elasticity |
| `enforceMonotonicEnergy` | false | true: fixed coefficients, no crawling, backtracking integrator |

**Event rule (shared)**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `useForceBasedTransitions` | true | Bell rule on; false = purely length-based T1 |
| `k_off0` | 1.0 | base Bell rate |
| `bellProbabilityThreshold` | 0.02 | threshold on `P_event` for T1 and division |

**Bulk T1 transitions**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `enableT1transitions` | true | T1 transitions on/off |
| `enhanceT1atWound` | true | use `T1_TOL_wound` for margin cells |
| `T1_TOL_bulk` / `T1_TOL_wound` | 0.05 / 0.15 | edge-length tolerances |
| `f_beta_T1` | 22.271433 | Bell force scale |
| `f_beta_length_slope_T1` | 1.0 | change of `f_beta_eff` with edge shortening |
| `T_eff_T1` | 0 | Metropolis temperature |
| `maxT1relaxsteps` | 30 | rest period of the vertices of an accepted T1 |

**Wound-margin handoffs and trims**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `enableWoundIntercalations` | true | margin events on/off |
| `L_intercalation_thresh` | 0.5 | absolute wound-contact length threshold |
| `L_intercalation_ratio_thresh` | 0.75 | contact length relative to the cell's mean edge length |
| `woundT1NewEdgeLength` | 0.01 | length of the new edge created by a handoff |
| `bellProbabilityThresholdT2` / `...T2Min` | 0.0125 / 0.009 | margin Bell threshold, upper / lower value |
| `f_beta_T2` | 0.494935 | Bell force scale |
| `f_beta_length_slope_T2` | 1.0 | change of `f_beta_eff` with contact shortening |
| `T_eff_T2` | 0 | Metropolis temperature |
| `minWoundBoundaryCells` | 5 | lower bound of the target margin population |
| `woundBoundaryCellEnergy` | 5.0 | energy per boundary cell in the acceptance test |

**Cell division (optional)**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `enableCellDivision` | false | division on/off |
| `allowWoundMarginDivision` | false | allow margin cells to divide |
| `Area_div_thresh` | 1.102315 | area threshold |
| `k_div0`, `f_beta_div`, `f_beta_div_slope` | 1.0, 0.250021, 1.0 | Bell parameters |
| `T_eff_div` | 3.512995 | Metropolis temperature |

**Wound sealing (optional)**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `enableWoundSealing` | false | sealing on/off |
| `L_seal_thresh` | 0.15 | vertex separation threshold |
| `k_seal0`, `f_beta_seal`, `f_beta_seal_slope` | 1.0, 2.2552, 1.0 | Bell parameters |
| `bellProbabilityThresholdSeal` | 0.02 | Bell threshold |
| `maxBoundaryCellsForSealing` | 5 | sealing starts at or below this margin size |
| `maxSealVertexDegree` | 5 | largest vertex degree allowed after a merge |
| `T_eff_seal` | 0 | Metropolis temperature |

**Plotting and other**

| Parameter | Default | Meaning |
| :--- | :--- | :--- |
| `nVisualisation` | 50 | plotting interval in the classical scripts |
| `plotForces`, `plotForces_Interval`, `plotForces_FigHandle`, `plotForces_Scale` | true, 20, 10, 2.0 | live force plot (`myMainClassical.m`) |
| `vel0`, `SelfPropellingCellIDs` | 0, 0 | self-propulsion (inactive) |
| `ApplyStretchX`, `StretchAtStep`, `StretchRatio`, `BoxIncompressibility` | 0, 100000, 1, 1 | uniaxial stretching (inactive) |
| `cellIDstoTrack` | 0 | cells highlighted in plots |
| `nT1`, `nWoundIntercalations`, `nWoundCornerTrims`, `nDivisions`, `nSeals` | 0 | event counters, updated during the run |

---

## 8. Output files

### 8.1 `runForceClosureStrategy.m`

| File | Content |
| :--- | :--- |
| `initial_state.png`, `final_state.png` | tissue at step 0 and at the last step, margin cells in magenta |
| `snapshots/wound_tNNNN.png` | tissue every 50 steps |
| `energy_vs_time.png` | mechanical energy vs model time |
| `wound_area_vs_time.png` | wound-loop area as % of its initial value |
| `wound_topology_diagnostics.png` / `.csv` | per-step loop area, area deficit, perimeter, boundary cells/edges, vertex degrees, smallest cell area and edge length, topology check |
| `wound_intercalation_events.csv` | per-step record of the margin-event evaluation (candidates, thresholds, accepted handoffs/trims, energy changes) |
| `wound_seal_events.csv` | per-step record of the sealing evaluation |
| `ka_contractility_every_50_steps.png` / `.csv` | margin area and contractility factors every 50 steps |
| `workspace.mat` | complete MATLAB workspace (parameters, final tissue, all logs) |
| `<run>.mp4` | video, when `WOUND_SAVE_VIDEO=1` |

Files added to the final results folder by `plot_contractility_annotations.m`: `contractility_vs_time.png/.fig`, `energy_vs_time_contractility_marked.png/.fig`, `contractility_vs_time_all_steps.csv`, `contractility_plot_notes.txt`.

### 8.2 `phaseDiagramValidation.m`

| File | Content |
| :--- | :--- |
| `phase_diagram_points.csv` | one row per relaxed point: set, Lbar, Gbar, p0, steps, status, mean shape index, edge-length variation, energy per cell, smallest area, phase |
| `phase_diagram.png` / `.fig` | simulated phases over the theoretical regions, and the shape index along the line |
| `relaxed_networks.png` | relaxed networks at two labelled points |

---

## 9. File reference

### 9.1 Simulation drivers

| File | Description |
| :--- | :--- |
| `runForceClosureStrategy.m` | Main wound-closure simulation used for the final results (Sections 2 and 5.1) |
| `runEnergyParameterTrial.m` | Function version of the main simulation with JSON overrides and a run tag; saves metrics and parameters as JSON |
| `continueEnergyTrial.m` | Continues a saved run for additional steps |
| `runMain.m` | Full closure run with video and closure metrics |
| `runForceBasedTransitions.m` | Edge-force-instrumented wound run with T1 force logs and Trepat-style profiles |
| `runEdgeForceAnalysis.m` | Wound run followed by edge-force decomposition plots |
| `myMainClassical.m` | Wound run with a live vertex-force plot |
| `myMainClassicalWoundOnly.m` | Wounded-tissue set-up and initial plot |
| `mainClassical.m` | Classical vertex model with T1 transitions, no wound |
| `phaseDiagramValidation.m` | Farhadifar phase-diagram validation with the model's own force, energy and update functions (Section 3.2) |

### 9.2 Parameter files

| File | Description |
| :--- | :--- |
| `SingleCellContractionParameters.m` | All wound-model parameters (Section 7) |
| `SelfPropelledParameters.m` | Parameters for self-propelled cell motion with the classical model |
| `StretchingParameters.m` | Parameters for uniaxial stretching with the classical model |

### 9.3 Mechanics

| File | Description |
| :--- | :--- |
| `getVertexForcesClassical.m` | Net vertex forces: area and perimeter terms (with margin multipliers), crawling, purse string, optional wound-area term |
| `getTissueEnergyClassical.m` | Energy used by the acceptance tests: cell energy, boundary-cell energy, purse-string energy, optional wound-area energy |
| `getAreaDerivative.m` | Analytical `dA/dr` of a cell with respect to one vertex |
| `getPerimeterDerivative.m` | Analytical `dP/dr` of a cell with respect to one vertex |
| `getEdgeForces.m` | Per-edge force data: perimeter tension, purse-string tension, area pressure, tangential and normal force vectors and their values per unit length |
| `updateVertexPositions.m` | Explicit Euler update with periodic wrapping |
| `updateVertexPositionsMonotonic.m` | Backtracking Euler update that keeps the energy non-increasing |

### 9.4 Topological events

| File | Description |
| :--- | :--- |
| `checkT1transitions.m` | Bulk T1 candidates, Bell threshold, trial swap and energy acceptance; logs the force at each accepted T1 |
| `getT1cellIDs.m` | The four cells around an edge selected for a T1 |
| `performT1swap.m` | Rewires connectivity and vertex-to-cell maps for a T1 |
| `checkWoundIntercalations.m` | Margin handoffs and corner trims: perimeter-paced target population, candidate classification, Bell threshold on compressive pressure, energy acceptance, event log |
| `performWoundMarginT1.m` | Executes a handoff (wound as the fourth cell) with a short perpendicular new edge |
| `trimWoundCorner.m` | Removes a corner vertex shared by two wound edges |
| `refreshWoundMargin.m` | Recomputes the free edges and margin cells after any topological change |
| `checkCellDivision.m` | Division candidates, Bell threshold on area pressure, energy acceptance |
| `performCellDivision.m` | Splits a cell through its centroid along its short axis |
| `getCellPrincipalAxis.m` | Long and short axes of a cell from its second-moment tensor |
| `insertVertexBetween.m` | Inserts a vertex between two cyclically adjacent vertices of a cell |
| `checkWoundSealing.m` | Sealing candidates (gap, Bell threshold on the closing force, energy acceptance) |
| `performWoundSeal.m` | Merges two margin vertices at their midpoint |

### 9.5 Wound geometry and diagnostics

| File | Description |
| :--- | :--- |
| `functions2D/createWound.m` | Removes cells inside a circular or elliptical region, using periodic barycenters and distances |
| `functions2D/findWoundEdgeCells.m` | Free edges (edges with one cell) and margin cells via a hash map of edges |
| `getOrderedWoundLoop.m` | Orders the free edges into a closed vertex loop |
| `getWoundDiagnostics.m` | Loop area, area deficit, perimeter, boundary counts, vertex degrees, smallest area and edge, topology check |
| `validateTissueTopology.m` | Checks finite coordinates, valid indices, non-degenerate cells and vertex-to-cell consistency |
| `getWoundArea.m` | Wound area from box-area conservation `Lx Ly - sum(A_c)` |
| `validateTrepatForces.m` | Radial and azimuthal edge-force profiles around the wound (single snapshot or time-pooled), with adaptive bins and error bars |

### 9.6 Cell geometry

| File | Description |
| :--- | :--- |
| `getCellAreas.m` / `getPolygonalCellArea.m` | Areas of all cells / of one cell |
| `getCellPerimeters.m` / `getPolygonalCellPerimeter.m` | Perimeters and edge lengths of all cells / of one cell |
| `getCellBarycenter.m` | Barycenter of a cell |
| `getTriangleArea.m` | Area of a triangle |
| `getCellDataWithinDomain.m` | Vertices and cells inside the domain |
| `getVertexToCellIDMap.m` | Vertex-to-cell map |
| `rearrangeAnticlockwise.m` | Orders a cell's vertices anticlockwise |
| `updateCellData.m` | Updates positions in `celldata` |

### 9.7 Mesh generation and periodicity

| File | Description |
| :--- | :--- |
| `genDataClassical.m` | Periodic Voronoi tissue from seeds placed with a minimum spacing; builds `celldata` |
| `checkInhibition.m` | Minimum-spacing test used when placing seeds |
| `checkSelfIntersect.m` | Vertex reordering for self-intersecting cells |
| `findMasters.m` | Maps periodic image vertices to their master vertices |
| `boundaryVertices.m` | Identifies boundary vertices through their images |
| `getTile.m` | 3 x 3 periodic tiling of seed points |
| `modifyVerticesForPeriodicity.m` | Unwraps a cell's vertices across the box edges |
| `updateImagePositions.m` | Moves image vertices with their masters |
| `updateMeshforPeriodicity.m` | Removes image vertices and rebuilds connectivity |

### 9.8 Visualisation and utilities

| File | Description |
| :--- | :--- |
| `plot3Dtissue.m` | Draws the cells as patches (margin and tracked cells highlighted) |
| `getCellDataforPlottingwithoutPeriJumps.m` | Prepares cell polygons for plotting across the periodic boundary |
| `plotEdgeForces.m` | Edge tension/pressure colour maps and force-vector plots |
| `PlotTissueEvolution.m` | Replays and plots a stored simulation |
| `PrintInformation.m` | Prints the step number and energy |
| `applyStretchX.m` | Applies uniaxial stretching along x |
| `insertelem.m` | Inserts an element into a vector |
| `cross3d.m` | Cross product of two 3-vectors |

### 9.9 `functions2D/`

Contains `createWound.m` and `findWoundEdgeCells.m`, together with the base versions of the shared helpers from the original classical model (`genDataClassical`, `getCellAreas`, `getTissueEnergyClassical`, `getVertexForcesClassical`, `performT1swap`, ...). When the scripts are run from `2D_VX_model/`, the files in `2D_VX_model/` take precedence over these base versions.

### 9.10 `results/`

| Folder | Content |
| :--- | :--- |
| `endsem_btp_results/` | Final wound-closure run (Section 3.1), including `plot_contractility_annotations.m` |
| `phase_diagram_validation_2026-11-19_185724/` | MATLAB phase-diagram validation (Section 3.2) |
| `phase_diagram_validation_numpy_port_2026-11-19/` | NumPy cross-check (Section 3.3): `vxport.py`, `sweep.py`, `plot_phase.py`, `sweep_results.csv`, `phase_diagram_validation.png`, `README.txt` |

### 9.11 Images

`30x30Circular1.png`, `30x30Circular2.png`: example wounded 30 x 30 tissues.

---

## 10. Main data structures

**`celldata`** (tissue state)

| Field | Content |
| :--- | :--- |
| `nCells`, `nMasterVertices` | number of cells and vertices |
| `r`, `r0`, `rinit` | current, previous and initial vertex positions (`nMasterVertices x 2`) |
| `connec` | cell array of anticlockwise vertex lists |
| `verttocell` | cell array: cells sharing each vertex |
| `A`, `P`, `Ainit`, `EdgeData` | areas, perimeters, initial areas, edge lengths per cell |
| `f` | vertex forces |
| `slave_r`, `StoMmap`, `idx*Boundary`, `boundaryNodes` | periodic image data |
| `woundEdges`, `cellIDtoContract` | free wound edges (`N x 2`) and margin cells |
| `woundArea0`, `woundPerimeter0` | initial wound-loop area and perimeter |
| `Kw_current`, `lambda_current`, `crawl_current`, `ka_wound_factor_current`, `contractility_wound_current` | recruitment state |

**`param`**: all parameters of Section 7, plus `Lx`, `Ly`, `woundType`, `woundCircular` (`center`, `radius`), `cellIDtoContract`, `initialWoundBoundaryCells` and the event counters.

**`edgeForceData`** (from `getEdgeForces.m`, one entry per unique edge): `v`, `cellIDs`, `length`, `edgeVec`, `midpoint`, `isWoundEdge`, `tension_peri`, `tension_purse`, `total_tension`, `pressure_area`, `force_tangent`, `force_normal`, `force_tangent_per_length`, `force_normal_per_length`.

---

## 11. References

1. Farhadifar, R., Roper, J.-C., Aigouy, B., Eaton, S. and Julicher, F. (2007). The influence of cell mechanics, cell-cell interactions, and proliferation on epithelial packing. *Current Biology*, 17, 2095-2104.
2. Brugues, A., Anon, E., Conte, V., Veldhuis, J. H., Gupta, M., Colombelli, J., Munoz, J. J., Brodland, G. W., Ladoux, B. and Trepat, X. (2014). Forces driving epithelial wound healing. *Nature Physics*, 10, 683-690.
3. Tetley, R. J., Staddon, M. F., Heller, D., Hoppe, A., Banerjee, S. and Mao, Y. (2019). Tissue fluidity promotes epithelial wound healing. *Nature Physics*, 15, 1195-1203.
4. Staddon, M. F., Bi, D., Tabatabai, A. P., Ajeti, V., Murrell, M. P. and Banerjee, S. (2018). Cooperation of dual modes of cell motility promotes epithelial stress relaxation to accelerate wound healing. *PLoS Computational Biology*, 14, e1006502.
5. Kaurin, D., Bal, P. K. and Arroyo, M. (2022). Peeling dynamics of fluid membranes bridged by molecular bonds: moving or breaking. *Journal of the Royal Society Interface*, 19, 20220183.
6. Bi, D., Yang, X., Marchetti, M. C. and Manning, M. L. (2016). Motility-driven glass and jamming transitions in biological tissues. *Physical Review X*, 6, 021011.
7. Trepat, X., Wasserman, M. R., Angelini, T. E., Millet, E., Weitz, D. A., Butler, J. P. and Fredberg, J. J. (2009). Physical forces during collective cell migration. *Nature Physics*, 5, 426-430.
