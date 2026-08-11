# Edge Force Calculation Strategy in the 2D Vertex Model

This document outlines the detailed mathematical strategy and algorithmic implementation used to calculate and visualize edge forces (tensions and pressures) in the `getEdgeForces.m` and `plotEdgeForces.m` scripts.

---

## 1. Theoretical Background: The Vertex Model Energy

In the classical 2D Vertex Model, the mechanical state of the epithelial tissue is governed by an energy functional for $N$ cells:

$$ E = \sum_{c=1}^{N} \left[ k_A (A_c - A_{c0})^2 + \frac{1}{2 r_{\text{stiff}}} (P_c - P_{c0})^2 \right] $$

Where for each cell $c$:
* $A_c$ is the current area, and $A_{c0}$ is the preferred area.
* $P_c$ is the current perimeter, and $P_{c0}$ is the preferred perimeter.
* $k_A$ is the area elasticity (resistance to volume changes).
* $r_{\text{stiff}}$ controls the perimeter stiffness (contractility vs adhesion).

The net force on any vertex $i$ is the negative gradient of this energy: $\mathbf{F}_i = -\nabla_{\mathbf{r}_i} E$. While standard simulations use this to move vertices, our goal for the edge force analysis is to attribute forces **to the edges themselves** to visualize how stress is distributed across the tissue.

---

## 2. Deriving Forces on an Edge

An edge $E_{ij}$ separates two cells (or belongs to one cell if at a boundary) and has length $L_{ij}$. The energy functional contributes two distinct types of mechanical stress across this edge:

### A. Perimeter Line Tension (Tangential Force)
The derivative of the perimeter term with respect to edge length $L_{ij}$ yields a scalar line tension $T_{\text{peri}}$ for cell $c$:
$$ T_{\text{peri}}^{(c)} = \frac{\partial E_{\text{peri}}}{\partial P_c} = \frac{2}{r_{\text{stiff}}} (P_c - P_{c0}) $$

Because an edge is shared by two adjacent cells (e.g., Cell 1 and Cell 2), the total line tension acting to contract (shorten) that edge is the sum of the contributions from both cells:
$$ T_{\text{peri}}^{\text{total}} = T_{\text{peri}}^{(1)} + T_{\text{peri}}^{(2)} $$

**Direction:** This force acts strictly **parallel** to the edge (tangentially).

### B. Area Pressure (Normal Force)
The area term acts like a 2D hydrostatic pressure inside the cell. The pressure inside cell $c$ is:
$$ \Pi^{(c)} = -\frac{\partial E_{\text{area}}}{\partial A_c} = -2 k_A (A_c - A_{c0}) $$

Note the negative sign: if a cell is compressed ($A_c < A_{c0}$), the pressure is positive (pushing outward). For an edge shared by Cell 1 and Cell 2, the net pressure difference across the edge exerts a normal force. 
In our formulation, we calculate the outward-pushing force contribution from each cell on the edge:
$$ F_{\text{pressure}}^{\text{total}} = \Pi^{(1)} L_{ij} + \Pi^{(2)} L_{ij} $$

**Direction:** This force acts strictly **perpendicular** (normal) to the edge.

### C. Wound Purse-String Contractility
If an edge lies on the wound margin, an active acto-myosin cable exerts an additional contractile line tension, denoted $\lambda_{\text{purse}}$.
$$ T_{\text{total}} = T_{\text{peri}}^{\text{total}} + \lambda_{\text{purse}} $$

---

## 3. Calculating Force Vectors (`getEdgeForces.m`)

To translate these scalar tensions and pressures into plottable 2D Cartesian vectors, we follow this algorithm:

1. **Edge Uniqueness:** We iterate over all cells and their vertices. To avoid processing shared edges twice, we use a Hash Map (`containers.Map`) keyed by the sorted vertex IDs (e.g., `'12_45'`). This allows us to accumulate the $T_{\text{peri}}$ and $\Pi$ contributions from all cells sharing that specific edge.
2. **Periodic Boundary Wrapping:** The edge vector $\mathbf{e} = \mathbf{r}_2 - \mathbf{r}_1$ is corrected for the periodic domain $[L_x, L_y]$:
   ```matlab
   edgeVec(1) = edgeVec(1) - Lx * round(edgeVec(1) / Lx);
   edgeVec(2) = edgeVec(2) - Ly * round(edgeVec(2) / Ly);
   ```
3. **Unit Vectors:** We define the unit tangent $\mathbf{\hat{t}}$ and unit normal $\mathbf{\hat{n}}$ (rotated 90° CCW from the tangent):
   $$ \mathbf{\hat{t}} = \frac{\mathbf{e}}{L_{ij}}, \quad \mathbf{\hat{n}} = [-\hat{t}_y, \hat{t}_x] $$
4. **Total Force Vectors:** 
   * **Tangential Force:** $\mathbf{F}_{\text{tan}} = T_{\text{total}} \cdot \mathbf{\hat{t}}$
   * **Normal Force:** $\mathbf{F}_{\text{norm}} = \text{Pressure} \cdot L_{ij} \cdot \mathbf{\hat{n}}$

### The New "Force Per Unit Length" Metric
While total normal force depends on the length of the edge ($\mathbf{F}_{\text{norm}} \propto L_{ij}$), the line tension does not. To compare force intensities fairly across edges of varying lengths (e.g., short edges near the wound vs. long edges far away), we compute the **Force per unit length**:
$$ \mathbf{f}_{\text{tan/length}} = \frac{\mathbf{F}_{\text{tan}}}{L_{ij}} $$
$$ \mathbf{f}_{\text{norm/length}} = \frac{\mathbf{F}_{\text{norm}}}{L_{ij}} = \text{Pressure} \cdot \mathbf{\hat{n}} $$

---

## 4. Visualization Strategy (`plotEdgeForces.m`)

Visualizing these forces effectively presents unique challenges:

1. **Midpoint Anchoring:** All quiver arrows originate from the edge midpoints (wrapped inside the periodic boundary).
2. **Auto-Scaling Prevention (The Alignment Fix):**
   * **The Problem:** MATLAB's `quiver(x, y, u, v, scale)` automatically normalizes and scales the `u,v` vectors based on the mean distance between `(x,y)` points. Because we plotted the tangential (green) and normal (blue) vectors in two separate `quiver` calls, MATLAB applied different auto-scale factors to them. This distorted the geometric angles, making the tangential arrows look "arbitrary" (not parallel to the edges).
   * **The Solution:** We explicitly disabled MATLAB's auto-scaling by passing `0` as the scale argument: `quiver(..., 0, ...)`.
3. **Manual Cross-Scaling:** 
   * Tangential forces per unit length ($T/L$) are numerically much smaller than normal forces ($\text{Pressure}$). 
   * To ensure both are visible simultaneously, we compute a specific `scaleRatio = max(\mathbf{f}_{\text{norm}}) / max(\mathbf{f}_{\text{tan}})`.
   * The tangential arrows are artificially boosted by this `scaleRatio` before plotting (only if they are significantly smaller). The legend reflects this multiplication (e.g., `Tangential/length (x15.2)`).
4. **Color Contrast:** We use High-contrast Green (`'g'`) for tangential forces and Blue (`'b'`) for normal forces so they remain distinctly visible against both the default yellow tissue cells and the red wound margin cells.
