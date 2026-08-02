# 2D Vertex Model (2D_VX_model)

This directory contains a MATLAB implementation of a **Classical 2D Vertex Model**, which is primarily used to simulate tissue mechanics, cellular evolution, and wound healing processes. The model relies on a pseudo-energy functional to calculate nodal forces and simulate vertex displacements over time. It also supports topological transitions (like T1 transitions) and periodic boundary conditions.

> **Note:** Many helper functions are located both in the root directory and the `functions2D/` subdirectory. Ensure that the MATLAB path is correctly set (e.g., `addpath('./functions2D')`).

---

## 1. Main Execution Scripts

These scripts are the entry points to the simulation. They initialize parameters, generate the starting tissue grid, and run the main time-stepping loop.

- **`mainClassical.m`**: The primary classical 2D vertex model simulating tissue mechanics with T1 transitions.
- **`myMainClassical.m`**: A modified version of the classical vertex model adapted for wound healing. It includes added features like force visualization during the main loop.
- **`myMainClassicalWoundOnly.m`**: Another entry point focused specifically on the wound healing scenario.

## 2. Model Parameters and Initialization

Scripts that define various specific physical configurations, forces, and initial generation of the cellular structure.

- **`SelfPropelledParameters.m`**: Defines parameters for self-propelled cell motion and related T1 transitions.
- **`SingleCellContractionParameters.m`**: Defines parameters to simulate the contraction of a single cell.
- **`StretchingParameters.m`**: Parameters used when a stretching force is applied to the tissue.
- **`genDataClassical.m`**: Generates the initial Voronoi tessellation to represent the cellular grid.

## 3. Core Physics, Energy, and Force Calculation

These scripts calculate the energy functional, derive the resulting forces on each vertex, and update the system's state based on dissipative dynamics (e.g., backward Euler method).

- **`getTissueEnergyClassical.m`**: Calculates the total pseudo-energy of the cellular collective based on cell areas and perimeters.
- **`getVertexForcesClassical.m`**: Calculates the forces acting on vertices by taking the gradient of the energy functional. Includes guard clauses to prevent simulation instability from collapsed cells.
- **`updateVertexPositions.m`**: Updates the current vertex positions using the calculated nodal forces, tissue viscosity, and time step.
- **`updateCellData.m`**: Updates cell data (areas, perimeters, etc.) based on the new positions of the nodes.

## 4. Geometric and Cellular Properties

Functions that calculate fundamental geometric properties for each cell in the tissue.

- **`getCellAreas.m` / `getPolygonalCellArea.m`**: Calculates the area of the cells (incorporating guards for invalid cells).
- **`getCellPerimeters.m` / `getPolygonalCellPerimeter.m`**: Calculates the perimeter of the cells and related edge data.
- **`getCellBarycenter.m`**: Finds the geometric center (barycenter) of a given cell.
- **`getCellDataWithinDomain.m`**: Extracts vertices and cell connectivity data that fall within a specified domain (for generating the relevant voronoi diagram).

## 5. Topological Transitions (T1 and Intercalations)

As tissue deforms, cell neighbors change. These scripts handle the topological T1 transitions and specific wound intercalations.

- **`checkT1transitions.m`**: Checks all edges in the tissue against a tolerance to see if a T1 transition (cell neighbor swap) should occur.
- **`getT1cellIDs.m`**: Identifies the four cells involved in a given T1 transition edge.
- **`performT1swap.m`**: Actually performs the neighbor swap by modifying nodal connectivity and mapping.
- **`checkWoundIntercalations.m`**: Handles cell detachments and specific intercalations at the wound margin (enhanced fluidity at the wound edge).

## 6. Wound Healing Specifics

Scripts that specifically deal with introducing a wound into the tissue and identifying its boundaries.

- **`functions2D/createWound.m`**: Removes a specific subset of cells to simulate a wound of a desired shape.
- **`functions2D/findWoundEdgeCells.m`**: An optimized hash-map-based function to find free edges in the tissue to identify the exact cells lining the wound boundary.

## 7. Periodic Boundary Conditions (PBC) & Tiling

To simulate an infinite or repeating tissue, periodic boundary conditions are applied. These functions map vertices that exit the bounding box back into it and handle "image" vertices.

- **`boundaryVertices.m`**: Identifies boundary vertices by checking against "slave" mappings.
- **`findMasters.m`**: Maps slave (image) vertices to their master (real) vertex IDs.
- **`getTile.m`**: Generates a 9x9 tiled representation of the tissue for periodic wrapping and visualization.
- **`modifyVerticesForPeriodicity.m`**: Updates vertex coordinates based on periodicity if they cross the domain bounds.
- **`updateImagePositions.m`**: Displaces image (slave) vertices exactly as their master vertices move.
- **`updateMeshforPeriodicity.m`**: Removes periodic images of vertices and reconstructs cell connectivity accordingly.

## 8. Plotting and Visualization

Tools used to visualize the tissue's state, forces, and evolution over time.

- **`PlotTissueEvolution.m`**: Main script for plotting the tissue as it evolves over time.
- **`plot3Dtissue.m`**: Plots the tissue using graphical patches. Includes fixes to avoid drawing invalid/collapsed cells.
- **`getCellDataforPlottingwithoutPeriJumps.m`**: Prepares cell coordinates for plotting by accounting for periodic boundary jumps (avoids drawing lines across the entire domain).
- **`PrintInformation.m`**: Helper for printing simulation progress or state information.

## 9. Utilities and Helpers

Miscellaneous helper functions for array manipulation and structural integrity.

- **`applyStretchX.m`**: Utility to apply a stretching force in the X direction.
- **`checkInhibition.m`**: Checks for node collisions based on an inhibition radius.
- **`checkSelfIntersect.m`**: Checks for and fixes internal element self-intersections (node swapping).
- **`insertelem.m`**: A generic utility to insert an element into a 1D vector at a given index.
- **`rearrangeAnticlockwise.m`**: Sorts a cell's vertices into an anti-clockwise order based on their angle relative to the barycenter.
