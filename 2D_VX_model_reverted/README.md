# 2D Vertex Model (Reverted Version) - `2D_VX_model_reverted`

This directory contains an older, "reverted" version of the **Classical 2D Vertex Model**. It is structurally very similar to the main `2D_VX_model` folder but contains specific variations in wound healing mechanics, force calculations, and execution scripts. 

> **Note:** Just like the main model, many helper functions are located both in the root directory and the `functions2D/` subdirectory.

---

## 1. Main Execution Scripts

- **`mainClassical.m`**: The primary classical 2D vertex model for basic tissue mechanics.
- **`my2MainClassical.m`**: A modified entry point for wound healing simulations. Features specific additions:
  - Visualizes forces inside the main time loop.
  - Calls a specific 3-cell intercalation function (`checkWoundIntercalations`).
  - Contains logic to save all output plots to a unique, timestamped folder to preserve parameters and runs.
- **`myMainClassicalWoundOnly.m`**: Another entry point focused specifically on the wound healing scenario.

## 2. Model Parameters and Initialization

- **`SelfPropelledParameters.m`**: Parameters for self-propelled cell motion and transitions.
- **`SingleCellContractionParameters.m`**: Parameters simulating contraction of individual cells.
- **`StretchingParameters.m`**: Parameters used when a stretching force is applied.
- **`genDataClassical.m`**: Generates the initial Voronoi tessellation to represent the cellular grid.

## 3. Core Physics, Energy, and Force Calculation

*Note: The force calculation here differs significantly from the newer version.*

- **`getTissueEnergyClassical.m`**: Calculates the total pseudo-energy of the cellular collective.
- **`getVertexForcesClassical.m`**: Calculates forces acting on vertices via energy gradients. **Important differences in this reverted version:**
  - *REMOVED:* The old whole-cell contraction approach (which modified the preferred perimeter $p_0$).
  - *ADDED:* A newer edge-specific "purse-string" force calculation, localizing contraction forces to specific wound edges.
- **`updateVertexPositions.m`**: Updates the current vertex positions using nodal forces.
- **`updateCellData.m`**: Updates cell geometry data based on new vertex positions.

## 4. Geometric and Cellular Properties

- **`getCellAreas.m` / `getPolygonalCellArea.m`**: Calculates cell areas.
- **`getCellPerimeters.m` / `getPolygonalCellPerimeter.m`**: Calculates cell perimeters and edge data.
- **`getCellBarycenter.m`**: Finds the geometric center of a cell.
- **`getCellDataWithinDomain.m`**: Extracts vertices and cells within a specified domain.

## 5. Topological Transitions & Specific Intercalations

- **`checkT1transitions.m`**: Checks for and performs T1 transitions. Uses a specific tolerance (`param.T1_TOL_wound`) for enhanced fluidity near wound margins.
- **`getT1cellIDs.m`**: Identifies cells involved in a T1 transition.
- **`performT1swap.m`**: Executes a neighbor swap. Takes a triggered T1 tolerance for vertex bounce mechanics.
- **`checkWoundIntercalations.m`**: **Important differences in this reverted version:**
  - Implements a 3-cell intercalation (topologically a T2 extrusion) at the wound margin using a "Vertex Merge & Orphan" method.
  - *Note:* This is "Version 1" of this logic and contains a known bug where the wound cell remains attached to the master vertex (`v1_master`).

## 6. Wound Healing Initialization (functions2D)

- **`functions2D/createWound.m`**: Removes specific cells based on periodic geometric centers to simulate a wound.
- **`functions2D/findWoundEdgeCells.m`**: A high-performance hash-map method to find free edges to map out the exact wound boundary.

## 7. Periodic Boundary Conditions (PBC) & Tiling

- **`boundaryVertices.m`**: Identifies boundary vertices.
- **`findMasters.m`**: Maps slave (image) vertices to their master IDs.
- **`getTile.m`**: Generates a 9x9 tiled representation of the tissue.
- **`modifyVerticesForPeriodicity.m`**: Updates coordinates when cells cross boundaries.
- **`updateImagePositions.m`**: Displaces image (slave) vertices exactly as their master vertices move.
- **`updateMeshforPeriodicity.m`**: Manages connectivity for periodic images.

## 8. Plotting and Visualization

- **`PlotTissueEvolution.m`**: Main script for plotting tissue evolution.
- **`plot3Dtissue.m`**: Plots the tissue using graphical patches.
- **`getCellDataforPlottingwithoutPeriJumps.m`**: Prepares cell coordinates for plotting by accounting for boundary jumps.
- **`PrintInformation.m`**: Helper for printing state information.

## 9. Utilities and Helpers

- **`applyStretchX.m`**: Applies stretching forces in the X direction.
- **`checkInhibition.m`**: Checks for node collisions.
- **`checkSelfIntersect.m`**: Checks for and fixes internal element self-intersections.
- **`insertelem.m`**: Helper to insert elements into a 1D vector.
- **`rearrangeAnticlockwise.m`**: Sorts a cell's vertices into an anti-clockwise order.
