# Force-Based Vertex Model of Epithelial Wound Closure

B.Tech. project, Department of Applied Mechanics, IIT Delhi.

A two-dimensional vertex model of epithelial wound closure in MATLAB. Polygonal cells move under area and perimeter elasticity, purse-string tension and lamellipodial crawling; neighbour exchanges and wound-margin rearrangements are selected from the local mechanical state and accepted through an energy criterion. The implemented mechanics are validated against the Farhadifar et al. (2007) soft/hexagonal phase diagram.

## Repository structure

| Folder | Content |
| :--- | :--- |
| `2D_VX_model/` | Current model, final results and phase-diagram validation. See [`2D_VX_model/README.md`](2D_VX_model/README.md) for the full description, how to reproduce the results and a reference for every file. |
| `2D_VX_model_reverted/` | Earlier development version of the model, kept for reference (see its own `README.md`). |

## Quick start

```matlab
cd('2D_VX_model')
runForceClosureStrategy      % main wound-closure simulation (results/ folder)
phaseDiagramValidation       % phase-diagram validation of the mechanics
```
