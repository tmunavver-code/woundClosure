# 2D Vertex Model for Epithelial Tissue Dynamics and Wound Closure

A biophysical MATLAB implementation of a **2D Vertex Model (VM)** designed to simulate epithelial tissue mechanics, cell shape dynamics, force transmission, and multi-mechanism epithelial wound closure. 

The model resolves explicit vertex and edge forces from a mechanical pseudo-energy functional, incorporates dissipative over-damped dynamics, handles periodic boundary conditions (PBC), and couples biological wound closure mechanisms (actomyosin purse-string contraction, lamellipodial cell crawling, dynamic margin stiffening, cell division, and force-gated topological transitions).

---

## Table of Contents
1. [Overview of the Model](#overview-of-the-model)
2. [Current Model State](#current-model-state)
3. [Biological & Biophysical Mechanisms](#biological--biophysical-mechanisms)
4. [File Directory and Intended Functions](#file-directory-and-intended-functions)
5. [How to Run the Code (Execution Guide)](#how-to-run-the-code-execution-guide)
   - [Prerequisites and Setup](#prerequisites-and-setup)
   - [Workflow 1: Complete Wound Closure Simulation (`runMain.m`)](#workflow-1-complete-wound-closure-simulation-runmainm)
   - [Workflow 2: Force-Based Closure Strategy Protocol (`runForceClosureStrategy.m`)](#workflow-2-force-based-closure-strategy-protocol-runforceclosurestrategym)
   - [Workflow 3: Force-Gated Topological Transitions & Stress Analysis (`runForceBasedTransitions.m`)](#workflow-3-force-gated-topological-transitions--stress-analysis-runforcebasedtransitionsm)
   - [Workflow 4: Edge Force Instrumentation & Trepat Monolayer Validation (`runEdgeForceAnalysis.m`)](#workflow-4-edge-force-instrumentation--trepat-monolayer-validation-runedgeforceanalysism)
   - [Workflow 5: Classical Vertex Model Baseline (`mainClassical.m`)](#workflow-5-classical-vertex-model-baseline-mainclassicalm)
   - [Workflow 6: Farhadifar Phase-Diagram Validation (`phaseDiagramValidation.m`)](#workflow-6-farhadifar-phase-diagram-validation-phasediagramvalidationm)
   - [Workflow 7: Parameterized / Continued Trials (`runEnergyParameterTrial.m`, `continueEnergyTrial.m`)](#workflow-7-parameterized--continued-trials-runenergyparametertrialm-continueenergytrialm)
6. [Key Model Parameters](#key-model-parameters)
7. [Simulation Outputs and Results](#simulation-outputs-and-results)

---

## Overview of the Model

In a 2D vertex model, an epithelial monolayer is represented as a planar polygonal network of vertices and shared edges. The mechanical state of the tissue is governed by a Farhadifar-style mechanical energy functional $E$:

$$E = \sum_{\alpha \in \text{cells}} \left[ k_a (A_\alpha - A_0)^2 + \frac{1}{r_{\text{stiff}}} (P_\alpha - p_0)^2 \right] + \sum_{\langle ij \rangle \in \text{edges}} \Lambda_{ij} L_{ij} + E_{\text{wound}}$$

Where:
- $A_\alpha$ and $P_\alpha$ are the area and perimeter of cell $\alpha$.
- $A_0$ and $p_0$ are the target (preferred) area and shape index ($p_0 = P_0 / \sqrt{A_0}$).
- $k_a$ represents bulk area compressibility modulus.
- $r_{\text{stiff}}$ controls perimeter elasticity and cortex contractility.
- $\Lambda_{ij}$ is the intercellular adhesion/line tension along edge $\langle ij \rangle$.
- $E_{\text{wound}}$ is an optional elastic-resistance term for the wound cavity (see [Current Model State](#current-model-state) — switched off by default in the current configuration).

The net conservative force acting on vertex $i$ is obtained from the gradient of the energy:

$$\mathbf{F}_i = -\nabla_{\mathbf{r}_i} E$$

Cell motion follows over-damped dissipative dynamics balancing conservative and active forces against substrate/intercellular viscous drag:

$$\eta \frac{d\mathbf{r}_i}{dt} = \mathbf{F}_i + \mathbf{F}_i^{\text{active}}$$

Vertex positions are iteratively updated using an explicit (forward) Euler scheme, $\mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) + \Delta t\, (\mathbf{F}_i + \mathbf{F}_i^{\text{active}}) / \eta$. A second, strictly energy-non-increasing integrator (`updateVertexPositionsMonotonic.m`) is also available: it backtracks the step size whenever a proposed move would raise the conservative energy, which is useful for runs immediately after a topological transition. All driver scripts seed MATLAB's random number generator (`rng('default')`) before building the initial mesh, so a given parameter set reproduces the same initial tessellation and simulation trajectory run to run.

---

## Current Model State

This section summarizes the current configuration of the model, as distinct from earlier versions of the code and the report.

- **Deterministic topological transitions.** T1 intercalations, wound-margin turnover (handoffs/trims), and wound sealing each still use a Bell's-law force-dependent rate to decide *eligibility*, but the final accept/reject step is now energy-gated deterministically: a transition is accepted only when it does not raise the total mechanical energy (the Metropolis uphill-acceptance temperatures `T_eff_T1`, `T_eff_T2`, and `T_eff_seal` are set to `0`). The stochastic Metropolis machinery itself remains in the code (and is still exercised by force-based cell division, which has its own non-zero `T_eff_div`), so a stochastic variant can be switched back on by raising the relevant `T_eff_*` parameter above zero.
- **Wound-hole restoring force removed from the active dynamics.** The elastic wound-hole resistance term ($K_w$, following Tetley et al., 2019) is still implemented in `getVertexForcesClassical.m` and `getTissueEnergyClassical.m`, but it is switched off by default (`param.enableWoundResistance = false`): it no longer exerts a vertex force or contributes to the energy. `K_w`'s exponential decay is retained purely as a timing clock that paces the recruitment of purse-string tension, crawling, and margin stiffening (via first-order relaxation toward decay-dependent targets, rate `k_recruit`). Wound closure is therefore driven entirely by the active purse-string and crawling forces, not by any mechanical "spring" pulling the hole back toward its original shape.
- **Farhadifar phase-diagram validation.** A new driver, `phaseDiagramValidation.m`, independently validates the vertex-model mechanics by reproducing the Farhadifar et al. (2007) soft/hexagonal network phase diagram directly from the repository's own `getVertexForcesClassical.m` / `getTissueEnergyClassical.m` functions — see [Workflow 6](#workflow-6-farhadifar-phase-diagram-validation-phasediagramvalidationm).
- **New diagnostics and sanity checks.** `getWoundDiagnostics.m` computes a per-step summary of wound geometry/topology (boundary cell/edge counts, wound perimeter, ordered loop area), and `validateTissueTopology.m` checks mesh connectivity each step (vertex-index bounds, finite coordinates, degenerate cells).
- **Parameterized and resumable trials.** `runEnergyParameterTrial.m` launches a fresh run from a JSON parameter-override file (`WOUND_CONFIG_JSON`) with a run tag (`WOUND_RUN_TAG`) and step count (`WOUND_NSTEPS`) read from environment variables, and `continueEnergyTrial.m` resumes a saved `workspace.mat` for additional steps under the same framework — see [Workflow 7](#workflow-7-parameterized--continued-trials-runenergyparametertrialm-continueenergytrialm).

---

## Biological & Biophysical Mechanisms

To capture realistic wound closure, the model couples physical force balances with biological responses:

1. **Actomyosin Purse-String Ring**: Cells along the wound margin recruit an active contractile actomyosin cable, applying additional line tension $\lambda$ along wound-facing edges.
2. **Active Lamellipodial Crawling**: Margin cells exert inward crawling traction forces directed toward the wound centroid.
3. **Kw Recruitment Clock**: A decaying variable $K_w(t) = K_{w0} e^{-t/\tau_{K_w}}$ (following Tetley et al., 2019) times the recruitment of purse-string tension, crawling, and margin stiffening toward their target values. In the current default configuration $K_w$ no longer exerts a wound-area elastic force or energy itself (`param.enableWoundResistance = false`) — that force term remains implemented and can be re-enabled — so it acts purely as a timing signal.
4. **Dynamic Margin Stiffening**: Margin cells undergo adaptive area compressibility and perimeter contractility stiffening ($k_a$ and $r_{\text{stiff}}$ ramp), recruited gradually toward their target multipliers (`ka_wound_factor`, `contractility_wound`) at the same pace as the Kw clock above.
5. **Force-Gated Topological Transitions (Bell's Law)**:
   - **T1 Transitions (Intercalations)**: Eligibility is driven by local edge tension $F$ via Bell's law rate:
     $$k_{\text{off}}(F) = k_{\text{off}0} \exp\left(\frac{F}{f_\beta}\right)$$
     In the current configuration the final accept/reject step is deterministic and energy-gated: a candidate transition is accepted only if it does not raise the total mechanical energy ($\Delta E \le 0$). The underlying Metropolis-style machinery (accepting an energy-increasing move with probability $\exp(-\Delta E / T_{\text{eff}})$) remains in the code for all transition types and is available by raising the corresponding `T_eff_*` parameter above zero.
   - **Margin Intercalations (Turnover)**: Cells bordering the margin exchange neighbors by treating the wound as a pseudo-cell (following Tetley et al., 2019), allowing cells to be expelled from the margin while remaining intact. The margin's target cell count is reduced gradually as the margin perimeter contracts, so handoffs track the shrinking wound rather than firing all at once.
   - **Wound Corner Trimming**: Non-convex corner vertices where both incident edges face the wound are trimmed to avoid topological deadlocks.
6. **Stress-Oriented Cell Division** *(available, off by default)*: Over-expanded cells divide when exceeding an area threshold, gated by Bell's law on area-elastic pressure. Division cleavage planes pass through the cell centroid along its short axis (orthogonal to the principal elongation/tensile axis). Enable with `param.enableCellDivision = true`.
7. **Wound Sealing** *(available, off by default)*: When opposing margins touch below a threshold distance, margin vertices coalesce, converting free boundary edges into internal edges and achieving complete scarless closure. Enable with `param.enableWoundSealing = true`.
8. **Farhadifar Phase-Diagram Validation**: Independently of the wound-closure mechanisms above, `phaseDiagramValidation.m` relaxes periodic hexagonal-lattice tissues across a sweep of normalized line tension/contractility and compares the resulting mean cell shape index against the theoretical Farhadifar et al. (2007) soft/hexagonal boundary ($p_0^* = \sqrt{8\sqrt{3}} \approx 3.7224$), using the same energy and force functions as the wound-closure simulations. See [Workflow 6](#workflow-6-farhadifar-phase-diagram-validation-phasediagramvalidationm).

---

## File Directory and Intended Functions

Below is the complete reference table of all files in `2D_VX_model/` and their specific functional roles.

| Filename | Category | Intended Function & Description |
| :--- | :--- | :--- |
| **`runMain.m`** | Driver Script | **Primary production simulation script.** Executes full 2000-timestep wound closure with purse-string, cell crawling, cell division, wound sealing, MP4 video generation, and comprehensive CSV metrics export to `final_results/`. |
| **`runForceClosureStrategy.m`** | Driver Script | *(Formerly `runSunday.m`)* 1000-timestep force-based closure experiment testing $K_w$ decay feedback, purse-string/crawling recruitment, and force-gated T1/T2 transitions. |
| **`runForceBasedTransitions.m`** | Driver Script | *(Formerly `runSaturday.m`)* Comprehensive force-instrumentation run validating Bell's law transitions, generating Trepat spatial profiles, force-vs-length scatter, and edge tension maps. |
| **`runEdgeForceAnalysis.m`** | Driver Script | Standalone analysis script that simulates tissue mechanics and produces 4-panel force decomposition plots and vector field maps for all cell edges. |
| **`phaseDiagramValidation.m`** | Driver / Validation | Reproduces the Farhadifar et al. (2007) soft/hexagonal phase diagram directly from the repository's energy/force functions: a force-vs-energy-gradient consistency check, a line sweep across $p_0$ at several $\bar\Gamma$ values, labelled example points, and relaxed-network figures. |
| **`runEnergyParameterTrial.m`** | Driver Script | Launches a fresh, reproducible run from a JSON parameter-override file and run tag, read from the `WOUND_CONFIG_JSON` / `WOUND_RUN_TAG` / `WOUND_NSTEPS` environment variables; supports the monotonic-energy integration mode. |
| **`continueEnergyTrial.m`** | Driver Script | Resumes a previously saved `workspace.mat` for additional timesteps under the same closure framework (used, e.g., for the `..._refined_eta2_extended` trials). |
| **`mainClassical.m`** | Driver Script | Baseline classical 2D vertex model driver (Sohan Kale, Adam Ouzeri, Nikhil Walani) simulating unperturbed tissue dynamics with standard T1 swaps. |
| **`myMainClassical.m`** | Driver Script | Classical vertex model variant adapted for wound scenarios with integrated real-time force visualization. |
| **`myMainClassicalWoundOnly.m`** | Driver Script | Minimal baseline script running classical wound relaxation without active closure mechanisms. |
| **`validateTrepatForces.m`** | Analysis / Validation | Validates spatial force profiles against Trepat et al. (2009/2014) monolayer stress microscopy profiles via radial and azimuthal force binning. |
| **`SingleCellContractionParameters.m`** | Parameter Config | Core parameter definitions: domain size ($L_x, L_y$), target shape index $p_0$, contractility, time stepping, wound geometry, Bell's law parameters, and division thresholds. |
| **`SelfPropelledParameters.m`** | Parameter Config | Configuration parameters for self-propelled Voronoi dynamics and active motility terms. |
| **`StretchingParameters.m`** | Parameter Config | Configuration parameters for external tissue stretching simulations along the X-axis. |
| **`genDataClassical.m`** | Mesh Generation | Generates the initial Voronoi cellular tessellation on a periodic domain with randomly seeded or regular cell centroids. |
| **`functions2D/createWound.m`** | Wound Geometry | Removes cells within a specified geometric profile (circular, elliptical, or rectangular) accounting for periodic boundaries to create the wound cavity. |
| **`functions2D/findWoundEdgeCells.m`** | Geometry / Mesh | High-performance $O(N)$ hash-map-based detector finding free boundary edges and the exact set of wound-margin cells. |
| **`getVertexForcesClassical.m`** | Physics Engine | Evaluates net conservative forces on all vertices via energy functional gradients ($-\nabla_{\mathbf{r}} E$). Includes stability guards for degenerate cells. |
| **`getEdgeForces.m`** | Physics Engine | Computes explicit tension and pressure acting along each unique edge: perimeter line tension, area pressure, purse-string tension, and force-per-length vectors. |
| **`getTissueEnergyClassical.m`** | Physics Engine | Calculates total mechanical energy $E$ of the tissue collective (area deviations, perimeter deviations, and active potentials). |
| **`getAreaDerivative.m`** | Physics Engine | Computes analytical derivative of cell area with respect to vertex coordinates $\partial A / \partial \mathbf{r}$. |
| **`getPerimeterDerivative.m`** | Physics Engine | Computes analytical derivative of cell perimeter with respect to vertex coordinates $\partial P / \partial \mathbf{r}$. |
| **`updateVertexPositions.m`** | Numerical Integrator | Updates vertex coordinates using explicit (forward) Euler integration based on nodal forces, viscosity $\eta$, and time step $\Delta t$. |
| **`updateVertexPositionsMonotonic.m`** | Numerical Integrator | Alternative integrator that backtracks the explicit-Euler step size until the conservative mechanical energy does not increase; selected via `param.enforceMonotonicEnergy`. |
| **`updateCellData.m`** | Numerical Integrator | Recomputes derived cell geometric properties (areas, perimeters) following vertex coordinate updates. |
| **`checkT1transitions.m`** | Topological Transitions | Detects short edges eligible for T1 neighbor swaps using Bell's law force gating and Metropolis acceptance. Differentiates bulk vs. wound margin tolerances. |
| **`getT1cellIDs.m`** | Topological Transitions | Identifies the four cells forming the quadrilateral tetrad around an edge undergoing a T1 swap. |
| **`performT1swap.m`** | Topological Transitions | Rewires mesh connectivity and vertex-to-cell maps to execute the topological T1 neighbor exchange. |
| **`checkWoundIntercalations.m`** | Topological Transitions | Evaluates 3-cell wound-margin intercalations (treating the wound as a pseudo-cell) to allow margin cell detachment and turnover. |
| **`performWoundMarginT1.m`** | Topological Transitions | Executes genuine T1-style neighbor swaps at the wound boundary, transferring vertices between margin neighbors so cells leave the margin intact. |
| **`trimWoundCorner.m`** | Topological Transitions | Removes sharp non-convex corner vertices where both incident edges are free wound boundaries, preventing margin turnover lockups. |
| **`checkWoundSealing.m`** | Topological Transitions | Evaluates whether opposing wound-margin vertices have met across the cavity gap to trigger force-gated adhesion and sealing. |
| **`performWoundSeal.m`** | Topological Transitions | Coalesces two colliding margin vertices into a single vertex, retiring boundary edges into standard internal edges to complete closure. |
| **`checkCellDivision.m`** | Cell Proliferation | Identifies enlarged cells exceeding area thresholds and triggers force-gated division via Bell's law stochastic pressure gating and energy acceptance. |
| **`performCellDivision.m`** | Cell Proliferation | Splits a dividing cell into two daughter cells along its short axis (orthogonal to principal elongation) through the centroid. |
| **`insertVertexBetween.m`** | Topology Utility | Inserts a newly generated vertex strictly between two adjacent vertices in a cyclic connectivity list, preserving edge orientation. |
| **`refreshWoundMargin.m`** | Topology Utility | Re-evaluates free edges and margin cell IDs directly from live connectivity after any topological event to guarantee consistency. |
| **`getCellAreas.m`** | Geometric Analysis | Computes areas for all cells in the tissue, returning zero for degenerated or eliminated cells. |
| **`getPolygonalCellArea.m`** | Geometric Analysis | Calculates the polygon area for a single cell using the shoelace formula with periodicity corrections. |
| **`getCellPerimeters.m`** | Geometric Analysis | Computes perimeters for all cells across the tissue mesh. |
| **`getPolygonalCellPerimeter.m`** | Geometric Analysis | Computes the boundary perimeter of a single cell polygon accounting for periodic edge wrapping. |
| **`getCellBarycenter.m`** | Geometric Analysis | Calculates the geometric center (barycenter) of a given cell polygon. |
| **`getCellPrincipalAxis.m`** | Geometric Analysis | Computes the principal elongation and short axes of a cell polygon using the second-moment shape tensor. |
| **`getWoundArea.m`** | Geometric Analysis | Computes total wound area robustly via box conservation ($A_{\text{wound}} = L_x L_y - \sum A_{\text{cell}}$), immune to margin pinching or loop self-intersections. |
| **`getOrderedWoundLoop.m`** | Geometric Analysis | Traverses unordered wound boundary edges to reconstruct a cyclic ordered vertex loop for margin geometric measurements. |
| **`getTriangleArea.m`** | Geometric Analysis | Elementary utility calculating the signed area of a 2D triangle from 3 vertices. |
| **`getCellDataWithinDomain.m`** | Geometric Analysis | Extracts vertices and connectivity of cells strictly inside the simulation domain bounds. |
| **`getWoundDiagnostics.m`** | Diagnostics | Per-step summary of wound geometry/topology: boundary cell/edge counts, wound perimeter, and ordered-loop area; used for logging and as the reference wound area/perimeter at run start. |
| **`validateTissueTopology.m`** | Diagnostics | Checks mesh connectivity and vertex-to-cell consistency each step (finite coordinates, in-range indices, non-degenerate cells). |
| **`boundaryVertices.m`** | Periodic Boundary Cond. | Identifies boundary vertices by matching periodic slave vertices to master coordinates. |
| **`findMasters.m`** | Periodic Boundary Cond. | Maps periodic slave (image) vertices along domain borders to their corresponding master vertex IDs. |
| **`modifyVerticesForPeriodicity.m`** | Periodic Boundary Cond. | Unwraps cell vertex coordinates crossing periodic boundary domain limits for accurate geometric evaluation. |
| **`updateImagePositions.m`** | Periodic Boundary Cond. | Enforces identical displacements on periodic image (slave) vertices as their master vertices move. |
| **`updateMeshforPeriodicity.m`** | Periodic Boundary Cond. | Strips periodic duplicate nodes and reconstructs the fundamental domain mesh. |
| **`getTile.m`** | Periodic Boundary Cond. | Generates a 3x3 (or 9-tile) representation of the periodic domain for boundary visualization. |
| **`plot3Dtissue.m`** | Visualization | Primary 2D graphical patch renderer plotting cells, wound-margin highlights, and periodic boundaries. |
| **`PlotTissueEvolution.m`** | Visualization | Generates multi-frame or time-lapse plots of tissue morphology evolution. |
| **`plotEdgeForces.m`** | Visualization | Renders color-coded scalar maps (line tension, area pressure) or quiver vector arrows on cell edges. |
| **`getCellDataforPlottingwithoutPeriJumps.m`** | Visualization | Eliminates spurious line artifacts across the periodic box by unwrapping polygons prior to rendering. |
| **`PrintInformation.m`** | Utilities | Formatted console logger printing step numbers, energies, and convergence info. |
| **`checkInhibition.m`** | Utilities | Spatial hash / collision check preventing vertices from being placed too close together upon initialization. |
| **`checkSelfIntersect.m`** | Utilities | Detects and resolves self-intersecting polygon vertices through barycentric sorting. |
| **`insertelem.m`** | Utilities | General 1D array insertion utility. |
| **`rearrangeAnticlockwise.m`** | Utilities | Orders cell vertex indices in counter-clockwise sequence around the cell barycenter. |
| **`cross3d.m`** | Utilities | Fast 3D cross-product calculation for 3-element vectors. |
| **`applyStretchX.m`** | Utilities | Applies a uniform uniaxial stretch displacement field to tissue vertices along the X axis. |
| **`final_results/`** | Results Folder | Output directory holding figures, simulation video (`wound_simulation.mp4`), CSV metrics, snapshots, and `workspace.mat`. |
| **`functions2D/`** | Function Library | Secondary folder containing core geometric and physics helper functions callable across simulations. |

---

## How to Run the Code (Execution Guide)

### Prerequisites and Setup
1. **MATLAB Version**: MATLAB R2020a or later is recommended.
2. **Toolboxes**: Basic MATLAB; Statistics and Machine Learning Toolbox (for Voronoi generation) and Image Processing Toolbox (optional, for frame capture).
3. **Working Directory**: Start MATLAB and set your current working directory to the `2D_VX_model/` folder:
   ```matlab
   cd('/path/to/woundClosure/2D_VX_model')
   ```
4. **Path Initialization**: Each script automatically adds `./functions2D` and `../common` to the MATLAB path:
   ```matlab
   addpath(genpath('./functions2D'))
   ```

---

### Workflow 1: Complete Wound Closure Simulation (`runMain.m`)
**Purpose**: The primary benchmark pipeline. Simulates a circular wound undergoing full closure over 2000 timesteps, combining purse-string tension, active cell crawling, margin contractility ramps, force-gated T1 transitions, cell divisions, and final wound sealing.

**Execution**:
In the MATLAB command prompt, type:
```matlab
runMain
```

**What it does**:
1. Initializes a periodic epithelial sheet of cells ($L_x = 15, L_y = 15$) and creates a central circular wound ($R = 2.0$).
2. Configures the $K_w$ recruitment clock, actomyosin purse-string, inward crawling traction, and margin contractility ramp (see [Current Model State](#current-model-state) for how these are currently coupled).
3. Steps through 2000 explicit-Euler integration steps with topological updates:
   - Bulk & margin T1 transitions (`checkT1transitions.m`)
   - Wound margin cell turnover (`checkWoundIntercalations.m`, `performWoundMarginT1.m`, `trimWoundCorner.m`)
   - Oriented cell division (`checkCellDivision.m`, `performCellDivision.m`)
   - Wound sealing (`checkWoundSealing.m`, `performWoundSeal.m`)
4. Renders and appends each timestep directly to an MP4 video (`wound_simulation.mp4`).
5. Saves snapshot figures every 50 timesteps into `final_results/snapshots/`.
6. Exports summary plots and numerical tables:
   - `wound_area_vs_time.png`
   - `margin_and_events_vs_time.png`
   - `energy_vs_time.png`
   - `ka_contractility_every_50_steps.png`
   - `closure_metrics.csv` (contains timestep, area %, margin cell count, edge count, cumulative T1, seals, and divisions)
   - `workspace.mat` (complete simulation workspace)

---

### Workflow 2: Force-Based Closure Strategy Protocol (`runForceClosureStrategy.m`)
**Purpose**: *(Formerly `runSunday.m`)* A 1000-timestep targeted experiment evaluating the dynamic feedback loop where the decay of the $K_w$ recruitment clock drives purse-string/crawling/contractility recruitment and force-gated topological transitions. This is the driver behind the project's designated results run. Run folders are named from the active mechanisms (e.g. a `with_resistance_` prefix appears only when `param.enableWoundResistance` is `true`), and an optional `WOUND_RUN_TAG` environment variable appends a custom tag to the folder name.

**Execution**:
```matlab
runForceClosureStrategy
```

**Outputs**:
Generates initial/final tissue states, energy evolution curves, wound area trajectory, contractility ramp CSV table, and saves the full workspace, under `results/TAHAFINAL/<descriptive-name>_<timestamp>/`.

---

### Workflow 3: Force-Gated Topological Transitions & Stress Analysis (`runForceBasedTransitions.m`)
**Purpose**: *(Formerly `runSaturday.m`)* Focuses on validating force-dependent junctional remodeling (Bell's law) against classical length-based transitions, generating stress maps and Trepat force comparisons.

**Execution**:
```matlab
runForceBasedTransitions
```

**Key Outputs**:
- Force-vs-edge-length scatter plots at T1 transition events.
- Azimuthal and radial spatial force profile comparing stress distribution to Trepat et al.
- Per-edge tension and pressure color maps.
- Step-by-step force vectors superimposed on the cell mesh.

---

### Workflow 4: Edge Force Instrumentation & Trepat Monolayer Validation (`runEdgeForceAnalysis.m`)
**Purpose**: Computes and decomposes edge forces across the tissue into normal area pressure, tangential line tension, and purse-string components.

**Execution**:
```matlab
runEdgeForceAnalysis
```
To validate against Trepat radial profiles using the workspace output:
```matlab
load('results/edge_forces/edge_force_workspace.mat');
validateTrepatForces(finalEdgeForces, param, param.Nsteps);
```

---

### Workflow 5: Classical Vertex Model Baseline (`mainClassical.m`)
**Purpose**: Runs the unperturbed classical 2D vertex model without wound closure modifications for baseline comparison.

**Execution**:
```matlab
mainClassical
```

---

### Workflow 6: Farhadifar Phase-Diagram Validation (`phaseDiagramValidation.m`)
**Purpose**: Independently validates that the tissue-mechanics implementation (`getVertexForcesClassical.m`, `getTissueEnergyClassical.m`) reproduces the known Farhadifar et al. (2007) soft/hexagonal phase diagram, by relaxing periodic hexagonal-lattice tissues across a sweep of normalized line tension $\bar\Lambda$ and contractility $\bar\Gamma$ and checking where the relaxed network's mean shape index departs from the hexagonal value $p_0^* = \sqrt{8\sqrt{3}} \approx 3.7224$.

**Execution**:
```matlab
phaseDiagramValidation
```
Set `runMode = 'quick'` (default, ~50 points) or `'full'` at the top of the script for a denser sweep.

**What it does**:
1. Runs a force-vs-energy-gradient consistency check (finite-difference vs. analytical force) on a perturbed hexagonal lattice and reports the relative error.
2. Relaxes a grid/line/point set of $(\bar\Lambda, \bar\Gamma)$ configurations to mechanical equilibrium and classifies each as hexagonal or soft from its relaxed shape index.
3. Writes `phase_diagram_points.csv` and renders the phase diagram (with the simulated points overlaid on the theoretical boundary) plus relaxed-network figures for selected points.

**Outputs** (saved to `results/phase_diagram_validation_<timestamp>/`):
- `phase_diagram_points.csv`
- `phase_diagram.png` / `.fig`
- `relaxed_networks.png`

---

### Workflow 7: Parameterized / Continued Trials (`runEnergyParameterTrial.m`, `continueEnergyTrial.m`)
**Purpose**: Run reproducible parameter sweeps and extend existing runs without editing `SingleCellContractionParameters.m` directly.

**Execution**:
```matlab
setenv('WOUND_CONFIG_JSON', '/path/to/overrides.json');
setenv('WOUND_RUN_TAG', 'my_trial');
setenv('WOUND_NSTEPS', '1500');          % optional
runEnergyParameterTrial('/path/to/overrides.json', 'my_trial')
```
`overrides.json` is a flat JSON object of `param` field names to override (the script errors on an unrecognized field name, to catch typos early). To extend a finished run instead:
```matlab
continueEnergyTrial('results/TAHAFINAL/<run_folder>/workspace.mat', 500)
```

---

## Key Model Parameters

Parameters are defined in [`SingleCellContractionParameters.m`](./SingleCellContractionParameters.m) and can be adjusted for sensitivity analysis, or overridden per-run via a JSON file as in [Workflow 7](#workflow-7-parameterized--continued-trials-runenergyparametertrialm-continueenergytrialm):

| Parameter | Default | Physical Meaning |
| :--- | :--- | :--- |
| `param.Lx`, `param.Ly` | `15.0, 15.0` | Periodic box dimensions (simulation domain size) |
| `param.p0` | `4.0` (fluid) / `2.0` (solid) | Target cell shape index ($P_0 / \sqrt{A_0}$) controlling tissue fluidity |
| `param.Nsteps` | `3000` (default); e.g. `1000` for `runForceClosureStrategy.m` | Number of numerical integration timesteps |
| `param.deltat` | `0.01` | Time increment per integration step |
| `param.enforceMonotonicEnergy` | `false` | When `true`, switches to the backtracking integrator (`updateVertexPositionsMonotonic.m`) that keeps the mechanical energy strictly non-increasing |
| `param.lambda_purse_string` | `4.0` | Actomyosin purse-string contractile tension, ramped up over `t_ramp_purse_string` (`50`) steps or by `Kw`-clock recruitment, depending on driver script |
| `param.crawl_force0` | `3.0` | Active lamellipodial crawling force per margin vertex; decays as the purse-string ramps up (`tau_decay_crawl` = `100` as a fallback time constant) |
| `param.enableWoundResistance` | `false` | Toggles whether $K_w$ exerts an actual elastic restoring force/energy on the wound opening. Off by default — closure is driven by the purse-string and crawling forces alone; `K_w` still runs as a recruitment-timing clock (see below) |
| `param.Kw0` | `2.0` | Initial value of the $K_w$ recruitment clock |
| `param.tau_Kw` | `2` | Characteristic decay timescale of the $K_w$ clock |
| `param.k_recruit` | `1.0` | Rate constant for the purse-string/crawling/margin-stiffening recruitment ODEs |
| `param.ka_wound_factor`, `param.contractility_wound` | `2.0`, `3.0` | Target multiplicative area-compressibility and contractility stiffening factors at the wound margin |
| `param.useForceBasedTransitions` | `true` | Enables Bell's law force-dependent T1/T2 transition eligibility |
| `param.f_beta_T1`, `param.f_beta_T2` | `22.27`, `0.49` | Bell's law characteristic bond force scales for T1 and wound-margin (T2-like) transitions |
| `param.T_eff_T1`, `param.T_eff_T2`, `param.T_eff_seal` | `0`, `0`, `0` | Metropolis uphill-acceptance temperature for T1, wound-margin turnover, and sealing. At `0`, acceptance is deterministic: a transition is accepted only if $\Delta E \le 0$ |
| `param.T_eff_div` | `3.51` | Metropolis uphill-acceptance temperature for cell division (non-zero; stochastic acceptance is available when `enableCellDivision` is turned on) |
| `param.enableCellDivision` | `false` | Master toggle for force-gated cell division (off by default) |
| `param.Area_div_thresh` | `1.10` | Relative cell area threshold triggering mitotic division |
| `param.enableWoundSealing` | `false` | Master toggle for wound-margin vertex coalescence/sealing (off by default) |
| `param.L_seal_thresh` | `0.15` | Opposing margin distance threshold triggering vertex coalescence (sealing) |

---

## Simulation Outputs and Results

When running the main simulation (`runMain.m`), all outputs are directed to the **`final_results/`** directory:

- **`wound_simulation.mp4`**: High-resolution video showing continuous tissue dynamics, margin cell tracking (highlighted in magenta), and wound closure progression.
- **`snapshots/`**: Contains numbered high-resolution PNG frames (`wound_t0050.png` through `wound_t2000.png`) sampled every 50 timesteps.
- **`closure_metrics.csv`**: Tabular data recording per-timestep wound area percentage, count of margin cells, number of free boundary edges, and cumulative event tallies (bulk T1 swaps, wound seals, cell divisions).
- **`wound_area_vs_time.png`**: Plot tracking wound area percentage over time from 100% ablation down to full closure.
- **`margin_and_events_vs_time.png`**: Dual-panel plot correlating the reduction in wound-margin cell count with cumulative topological events.
- **`ka_contractility_every_50_steps.csv` / `.png`**: Records the dynamic ramp of area modulus ($k_a$) and contractility ($1/r_{\text{stiff}}$) along the wound margin.
- **`energy_vs_time.png`**: Total mechanical pseudo-energy trajectory verifying numerical dissipation and stability.
- **`workspace.mat`**: Complete MATLAB binary workspace holding all variables, parameters, and time series for post-hoc analysis.

`runForceClosureStrategy.m`, `runEnergyParameterTrial.m`, and `continueEnergyTrial.m` instead save each run to its own tagged, timestamped folder under **`results/TAHAFINAL/`**, with the same style of outputs (plus an optional every-timestep video via `WOUND_SAVE_VIDEO=1`), so parameter sweeps and extended trials don't overwrite one another.

`phaseDiagramValidation.m` saves independently to **`results/phase_diagram_validation_<timestamp>/`**: `phase_diagram_points.csv`, `phase_diagram.png` / `.fig`, and `relaxed_networks.png` (see [Workflow 6](#workflow-6-farhadifar-phase-diagram-validation-phasediagramvalidationm)).

---

### Reference Publications
1. **Farhadifar, R., et al.** (2007). *The influence of cell mechanics, cell-cell interactions, and cell division on epithelial packing*. Current Biology, 17(24), 2095-2104.
2. **Trepat, X., et al.** (2009). *Physical forces during collective cell migration*. Nature Physics, 5(6), 426-430.
3. **Tetley, R. J., et al.** (2019). *Tissue fluidification by cell shape changes during embryonic repair*. Nature Physics, 15(11), 1195-1203.
