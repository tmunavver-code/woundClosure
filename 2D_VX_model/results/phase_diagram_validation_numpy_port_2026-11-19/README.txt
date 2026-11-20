Phase-diagram validation of the 2D_VX_model mechanics (Farhadifar et al. 2007)
=============================================================================

Mapping (this model -> Farhadifar):
  Gbar = 1/(ka*rstiff)      Lbar = -2*p0/(ka*rstiff) = -2*p0*Gbar
  soft/hexagonal boundary  = the ray p0 = sqrt(8*sqrt(3)) = 3.7224

Method: periodic hexagonal lattice of unit-area cells (fixed box, as in all
repo simulations), perturbed by 0.02, relaxed with the repo's force law
(explicit Euler, dt = 0.01, eta = 1, ka = 1), no wound, no T1s.
Classified by mean shape index q = <P/sqrt(A)>: hexagonal if q - 3.7224 <= 1e-3.

Results (sweep_results.csv, 285 runs):
  - forces = -grad(energy) of the repo's own functions to ~1e-10 (relative)
  - along the line, all six Gbar rows (0.02-0.17): hexagonal at p0 = 3.720,
    soft at p0 = 3.725 -> boundary within +-0.001 in Lbar of theory
  - coarse grid (189 points): 0 disagreements with theory
  - labelled points I, II, V hexagonal; III, IV soft
  - hatched "collapse" region needs a pressure-free box; not tested
    (in the fixed box the hexagonal state stays stable to Lbar = 2.7)

These numbers come from a NumPy port (vxport.py) of getPolygonalCellArea,
getAreaDerivative, getPerimeterDerivative, getVertexForcesClassical,
getTissueEnergyClassical and updateVertexPositions, because MATLAB was not
available where they were run. To confirm on the MATLAB code itself, run
2D_VX_model/phaseDiagramValidation.m (it calls those functions directly).

Files: phase_diagram_validation.png, sweep_results.csv,
       vxport.py (port), sweep.py (runs), plot_phase.py (figure)
