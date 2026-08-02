# BTP_Report

## Project Structure
This repository contains two main directories for the 2D vertex models simulating tissue mechanics and wound healing: `2D_VX_model` and `2D_VX_model_reverted`. 

### Differences between `2D_VX_model` and `2D_VX_model_reverted`

While both folders contain similar codebase foundations for vertex model simulation, they represent different stages of development with critical differences in physics implementation and stability:

1. **Topological Transitions (`checkWoundIntercalations.m`)**
   - **`2D_VX_model`**: Uses "VERSION 2.1 (Corrected)" which properly detaches a cell from the wound margin during a 3-cell intercalation (extrusion).
   - **`2D_VX_model_reverted`**: Uses "VERSION 1", which contains a known bug where extruded cells fail to fully detach and remain tethered to a master vertex.

2. **Simulation Stability (`getVertexForcesClassical.m`)**
   - **`2D_VX_model`**: Includes a "guard clause" that skips force calculations for degenerated/collapsed cells (cells with fewer than 3 vertices). This prevents the physics engine from dividing by zero and breaking the simulation.
   - **`2D_VX_model_reverted`**: Missing this stability fix, making it prone to crashing when cells collapse.

3. **Rendering & Plotting (`plot3Dtissue.m`)**
   - **`2D_VX_model`**: Includes a fix for a "white hole" bug, successfully ignoring the drawing step for invalid polygons (fewer than 3 vertices).
   - **`2D_VX_model_reverted`**: Attempts to draw all polygons indiscriminately, which can cause rendering glitches when invalid cells exist.

4. **Execution & Saving Logic**
   - **`2D_VX_model`**: Uses `myMainClassical.m` as the main script. It is more straightforward and was historically used for direct visualization. 
   - **`2D_VX_model_reverted`**: Uses `my2MainClassical.m` which features a more sophisticated, automated saving system. It dynamically generates timestamped folders incorporating the physical parameters (e.g., `2026-08-02_143948_fluid_Circ_r1.5_...`) and saves plots invisibly in the background.

**Conclusion**: Use `2D_VX_model` for physically accurate and stable wound healing simulations. The `2D_VX_model_reverted` version has a better automated saving structure but relies on broken physics logic for cell extrusions.