Phase-diagram validation of the 2D_VX_model mechanics (Farhadifar et al. 2007)
Independent NumPy reimplementation
=============================================================================

Mapping (this model -> Farhadifar):
  Gbar = 1/(ka*rstiff)      Lbar = -2*p0/(ka*rstiff) = -2*p0*Gbar
  soft/hexagonal boundary  = the ray p0 = sqrt(8*sqrt(3)) = 3.7224

Method: periodic hexagonal lattice of unit-area cells (fixed box, as in the
wound simulations), perturbed by 0.02, relaxed with the model's force law
(explicit Euler, dt = 0.01, eta = 1, ka = 1), no wound, no T1 transitions.
Classified by mean shape index q = <P/sqrt(A)>: hexagonal if q - 3.7224 <= 1e-3.

Results (sweep_results.csv, 285 relaxations):
  - vertex forces equal the negative energy gradient to ~1e-9 (relative)
  - along the line, all six Gbar rows (0.02-0.17): hexagonal at p0 = 3.720,
    soft at p0 = 3.725, bracketing the theoretical boundary
  - coarse grid (189 points): every point matches the theoretical phase
  - labelled points I, II, V hexagonal; III, IV soft
  - in the fixed box the hexagonal state stays stable up to Lbar = 2.7

vxport.py reimplements getPolygonalCellArea, getAreaDerivative,
getPerimeterDerivative, getVertexForcesClassical, getTissueEnergyClassical and
updateVertexPositions in NumPy. The MATLAB counterpart of this test is
2D_VX_model/phaseDiagramValidation.m, which calls the MATLAB functions directly.

Files:
  vxport.py                     NumPy implementation of the bulk mechanics
  sweep.py                      runs all relaxations, writes sweep_results.csv
                                and the example networks (cfg_*.npz)
  plot_phase.py                 draws phase_diagram_validation.png
  sweep_results.csv             one row per relaxation
  phase_diagram_validation.png  phase diagram and line-crossing plot

To regenerate (Python 3 with numpy, pandas, matplotlib):
  python3 sweep.py
  python3 plot_phase.py
