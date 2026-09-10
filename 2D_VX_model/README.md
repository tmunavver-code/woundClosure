# 2D Vertex Model for Epithelial Tissue Dynamics and Wound Closure

A biophysical MATLAB implementation of a **2D Vertex Model (VM)** designed to simulate epithelial tissue mechanics, cell shape dynamics, force transmission, and multi-mechanism epithelial wound closure. 

The model resolves explicit vertex and edge forces from a mechanical pseudo-energy functional, incorporates dissipative over-damped dynamics, handles periodic boundary conditions (PBC), and couples biological wound closure mechanisms (actomyosin purse-string contraction, lamellipodial cell crawling, dynamic margin stiffening, cell division, and force-gated topological transitions).

---

## Table of Contents
1. [Overview of the Model](#overview-of-the-model)
2. [Biological & Biophysical Mechanisms](#biological--biophysical-mechanisms)
3. [File Directory and Intended Functions](#file-directory-and-intended-functions)
4. [How to Run the Code (Execution Guide)](#how-to-run-the-code-execution-guide)
   - [Prerequisites and Setup](#prerequisites-and-setup)
   - [Workflow 1: Complete Wound Closure Simulation (`runMain.m`)](#workflow-1-complete-wound-closure-simulation-runmainm)
   - [Workflow 2: Force-Based Closure Strategy Protocol (`runForceClosureStrategy.m`)](#workflow-2-force-based-closure-strategy-protocol-runforceclosurestrategym)
   - [Workflow 3: Force-Gated Topological Transitions & Stress Analysis (`runForceBasedTransitions.m`)](#workflow-3-force-gated-topological-transitions--stress-analysis-runforcebasedtransitionsm)
   - [Workflow 4: Edge Force Instrumentation & Trepat Monolayer Validation (`runEdgeForceAnalysis.m`)](#workflow-4-edge-force-instrumentation--trepat-monolayer-validation-runedgeforceanalysism)
   - [Workflow 5: Classical Vertex Model Baseline (`mainClassical.m`)](#workflow-5-classical-vertex-model-baseline-mainclassicalm)
5. [Key Model Parameters](#key-model-parameters)
6. [Simulation Outputs and Results](#simulation-outputs-and-results)

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
- $E_{\text{wound}}$ accounts for the elastic resistance of the wound cavity.

The net conservative force acting on vertex $i$ is obtained from the gradient of the energy:

$$\mathbf{F}_i = -\nabla_{\mathbf{r}_i} E$$

Cell motion follows over-damped dissipative dynamics balancing conservative and active forces against substrate/intercellular viscous drag:

$$\eta \frac{d\mathbf{r}_i}{dt} = \mathbf{F}_i + \mathbf{F}_i^{\text{active}}$$

Using a backward Euler time integration scheme, vertex positions are iteratively updated to advance the tissue state.

---

## Biological & Biophysical Mechanisms

To capture realistic wound closure, the model couples physical force balances with biological responses:

1. **Actomyosin Purse-String Ring**: Cells along the wound margin recruit an active contractile actomyosin cable, applying additional line tension $\lambda$ along wound-facing edges.
2. **Active Lamellipodial Crawling**: Margin cells exert inward crawling traction forces directed toward the wound centroid.
3. **Decaying Wound Hole Resistance ($K_w$)**: Freshly ablated wounds exhibit initial matrix/fluid resistance that decays over characteristic time $\tau_{K_w}$, dynamically modulating contractility recruitment.
4. **Dynamic Margin Stiffening**: Margin cells undergo adaptive area compressibility and perimeter contractility stiffening ($k_a$ and $r_{\text{stiff}}$ ramp) gated by wound resistance feedback.
5. **Force-Gated Topological Transitions (Bell's Law)**:
   - **T1 Transitions (Intercalations)**: Junctional remodeling is driven stochastically by local edge tension $F$ via Bell's law rate:
     $$k_{\text{off}}(F) = k_{\text{off}0} \exp\left(\frac{F}{f_\beta}\right)$$
     with a Metropolis-style energy acceptance criterion ($\Delta E \le 0$ or accepted with probability $\exp(-\Delta E / T_{\text{eff}})$).
   - **Margin Intercalations (Turnover)**: Cells bordering the margin exchange neighbors by treating the wound as a pseudo-cell (following Tetley et al., 2019), allowing cells to be expelled from the margin while remaining intact.
   - **Wound Corner Trimming**: Non-convex corner vertices where both incident edges face the wound are trimmed to avoid topological deadlocks.
6. **Stress-Oriented Cell Division**: Over-expanded cells divide when exceeding an area threshold. Division cleavage planes pass through the cell centroid along its short axis (orthogonal to the principal elongation/tensile axis).
7. **Wound Sealing**: When opposing margins touch below a threshold distance, margin vertices coalesce, converting free boundary edges into internal edges and achieving complete scarless closure.

---

## File Directory and Intended Functions

Below is the complete reference table of all files in `2D_VX_model/` and their specific functional roles.

| Filename | Category | Intended Function & Description |
| :--- | :--- | :--- |
| **`runMain.m`** | Driver Script | **Primary production simulation script.** Executes full 2000-timestep wound closure with purse-string, cell crawling, cell division, wound sealing, MP4 video generation, and comprehensive CSV metrics export to `final_results/`. |
| **`runForceClosureStrategy.m`** | Driver Script | *(Formerly `runSunday.m`)* 1000-timestep force-based closure experiment testing $K_w$ decay feedback, purse-string/crawling recruitment, and force-gated T1/T2 transitions. |
| **`runForceBasedTransitions.m`** | Driver Script | *(Formerly `runSaturday.m`)* Comprehensive force-instrumentation run validating Bell's law transitions, generating Trepat spatial profiles, force-vs-length scatter, and edge tension maps. |
| **`runEdgeForceAnalysis.m`** | Driver Script | Standalone analysis script that simulates tissue mechanics and produces 4-panel force decomposition plots and vector field maps for all cell edges. |
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
| **`updateVertexPositions.m`** | Numerical Integrator | Updates vertex coordinates using backward Euler integration based on nodal forces, viscosity $\eta$, and time step $\Delta t$. |
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
2. Configures decaying hole resistance ($K_w$), actomyosin purse-string, inward crawling traction, and margin contractility ramp.
3. Steps through 2000 backward-Euler integration steps with topological updates:
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
**Purpose**: *(Formerly `runSunday.m`)* A 1000-timestep targeted experiment evaluating the dynamic feedback loop where the decay of wound hole resistance ($K_w$) drives contractility recruitment and force-gated topological transitions.

**Execution**:
```matlab
runForceClosureStrategy
```

**Outputs**:
Generates initial/final tissue states, energy evolution curves, wound area trajectory, contractility ramp CSV table, and saves the full workspace.

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

## Key Model Parameters

Parameters are defined in [`SingleCellContractionParameters.m`](file:///Users/tmunavver/Downloads/woundClosure/2D_VX_model/SingleCellContractionParameters.m) and can be adjusted for sensitivity analysis:

| Parameter | Default | Physical Meaning |
| :--- | :--- | :--- |
| `param.Lx`, `param.Ly` | `15.0, 15.0` | Periodic box dimensions (simulation domain size) |
| `param.p0` | `3.72` (fluid) / `3.81` (solid) | Target cell shape index ($P_0 / \sqrt{A_0}$) controlling tissue fluidity |
| `param.Nsteps` | `1000` / `2000` | Number of numerical integration timesteps |
| `param.dt` | `0.01` | Time increment per integration step |
| `param.lambda_wound_max` | `0.40` | Maximum actomyosin purse-string contractile tension |
| `param.crawl_force_max` | `0.08` | Maximum active lamellipodial crawling force per margin vertex |
| `param.Kw0` | `2.0` | Initial elastic resistance of the unclosed wound cavity |
| `param.tau_Kw` | `300` | Characteristic timescale for wound resistance decay |
| `param.ka_wound_ramp_max` | `3.0` | Multiplicative area-compressibility stiffening factor at margin |
| `param.useForceBasedTransitions` | `true` | Enables Bell's law force-dependent T1 transitions |
| `param.f_beta` | `0.50` | Bell's law characteristic bond force scale |
| `param.Area_div_thresh` | `1.30` | Relative cell area threshold triggering mitotic division |
| `param.L_seal_thresh` | `0.40` | Opposing margin distance threshold triggering vertex coalescence (sealing) |

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

---

### Reference Publications
1. **Farhadifar, R., et al.** (2007). *The influence of cell mechanics, cell-cell interactions, and cell division on epithelial packing*. Current Biology, 17(24), 2095-2104.
2. **Trepat, X., et al.** (2009). *Physical forces during collective cell migration*. Nature Physics, 5(6), 426-430.
3. **Tetley, R. J., et al.** (2019). *Tissue fluidification by cell shape changes during embryonic repair*. Nature Physics, 15(11), 1195-1203.
