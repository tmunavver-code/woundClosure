# Research Papers Summary

This file summarizes the three papers used as the theoretical basis for the
force-based wound-closure strategy implemented in this directory. Each
section covers what the paper actually says; the final section maps that
content onto the specific files and parameters in this codebase.

---

## 1. Trepat et al. 2014 — "Forces driving epithelial wound healing"

**Citation:** Brugués A, Anon E, Conte V, Veldhuis JH, Gupta M, Colombelli J,
Muñoz JJ, Brodland GW, Ladoux B, Trepat X. *Nature Physics* 10, 683–690
(2014). `Xavi_NP_2014.pdf`

### What they did
Ablated a ~5,000 μm² cluster of ~20 cells in a confluent MDCK monolayer
(LifeAct-GFP, on collagen-coated polyacrylamide gel with fluorescent bead
markers) and measured, simultaneously, cell velocities (PIV), substrate
traction forces (traction-force microscopy and flexible micropillar
arrays), and actomyosin/focal-adhesion structure (immunofluorescence for
actin, phospho-myosin, paxillin, talin) throughout closure.

### Key findings

- **Two temporally distinct force mechanisms**, not one:
  - **OPTL (outward-pointing traction layer):** in the first row of cells at
    the leading edge, early in closure, cells extend lamellipodia/filopodia
    into the wound and generate tractions on the substrate that point
    *away* from the wound (i.e., outward) — the classic signature of
    lamellipodial crawling.
  - **IPTL (inward-pointing traction layer):** appearing shortly behind the
    leading edge and growing with time, a second traction layer pointing
    *toward* the wound. This was novel — not previously reported.
- The **supracellular actomyosin ring (purse-string)** forms with a ~15 min
  delay after ablation and remains present throughout, coexisting with
  ongoing lamellipodial protrusion (the two mechanisms are not mutually
  exclusive, contrary to older either/or framing).
- **The IPTL originates from the ring, not the leading-edge cells' own
  crawling.** Evidence: focal adhesions accumulate underneath the ring
  (basally) as closure proceeds and are oriented *tangentially* to the
  wound edge, not radially — this is what actually explains the tangential
  component of IPTL traction, since a ring under uniform tension alone
  would only produce radial (inward) substrate displacement, not tangential
  traction. Ablating the ring itself causes an immediate, local drop in
  IPTL tractions, confirming the ring is the source.
- **Calcium chelation (EGTA)** prevents ring assembly (junctions weakened
  but not destroyed) — cells still close the wound at a similar rate via
  crawling alone, but the IPTL never appears and directionality of motion
  is lost. This decouples the two mechanisms experimentally.
- **Force pattern is insensitive to substrate stiffness** (tested down to
  3 kPa, threefold softer) even though single-cell migration on soft
  substrates is known to slow down — collective wound closure is more
  robust than single-cell motility to this variable.
- Built a matching 2D finite-element continuum model (viscous cytoplasm +
  incompressible) combining lamellipodial edge forces and actomyosin-ring
  tension; reproduces the traction patterns only when *both* mechanisms and
  the ring's *tangential* force component are included.

### Relevance to this codebase
Motivates the **wound-directed crawling force** (`param.crawl_force0`,
`param.tau_decay_crawl` in `getVertexForcesClassical.m`) as an early, active,
independent driver on wound-margin cells, separate from and temporally
prior to the purse-string ring's dominance — the OPTL→IPTL handoff. Note the
important caveat already recorded elsewhere in this repo: **this vertex
model has no substrate**, so implementing the crawling force as a direct
pull toward the wound centroid is a pattern-level analogy to Trepat's
mechanism, not a literal reproduction of traction-force magnitudes or the
focal-adhesion-mediated tangential-force mechanism they identify for the
IPTL.

---

## 2. Tetley et al. 2019 — "Tissue fluidity promotes epithelial wound healing"

**Citation:** Tetley RJ, Staddon MF, Heller D, Hoppe A, Banerjee S, Mao Y.
*Nature Physics* 15, 1195–1203 (2019). Supplementary Information used here:
`SI_Tissue_Fluidity.pdf` (Methods, Supp. Figs. 1–10, Supp. Table 1, Video
legends — the main-text figures/results were not separately provided, but
the Methods section reproduces the full model and most quantitative claims
needed here).

### Experimental system
*Drosophila* wing imaginal discs, laser-ablated to create wounds, imaged
live (confocal, 3 min intervals) and tracked/segmented with EpiTools.
Genotypes compared: wild-type (WT), *Rok*-RNAi (Rho-kinase knockdown,
reduces actomyosin contractility → more fluid), *Mbs*-RNAi (myosin
phosphatase knockdown, increases contractility → more solid/jammed), and
ubi-Ecad-GFP (uniform E-cadherin labeling, used for the row-resolved
intercalation analysis). Also embryos (Ecad-tdTomato) for a cross-species
intercalation comparison.

### The vertex model (their formulation — important: different from this
repo's parametrization, see the mapping note in section 4 below)

```
E = Σ_α (1/2) K (A_α - A_α^0)^2  +  Σ_α (1/2) Γ P_α^2  +  Σ_⟨i,j⟩ Λ_ij L_ij
```

- Area elasticity (`K`), **contractility** (`Γ`, quadratic in perimeter —
  no separate target-perimeter term), and an explicit **line tension**
  `Λ_ij` on every edge, evolved dynamically.
- Overdamped dynamics: `μ dx_i/dt = F_i = -∂E/∂x_i`.
- Shape index `p0 = -Λ0 / (2Γ)` is the fluid/solid jamming control
  parameter (from Bi et al. 2015), with the transition at **`p0* ≈ 3.81`**.
  Their own shear/viscosity sweep (Supp. Fig. 5e-f) shows both quantities
  decreasing toward ~0 as `p0` rises past 3.8.

**Two sources of activity keep the tissue out of equilibrium:**
1. **Cell division** — resting/growing/dividing cycle; preferred area
   doubles over 30 min in mitosis, then splits into two cells along the
   short axis at the lower-energy configuration. Mean division rate fit to
   experimental measurements (~2.25%/hr WT).
2. **Line tension fluctuations** — an Ornstein-Uhlenbeck process on every
   edge, bulk included: `dΛ_ij/dt = -(1/τ_m)(Λ_ij - Λ0) + ξ_ij(t)`, with
   white noise `⟨ξ_ij(t)ξ_kl(t')⟩ = (2σ_m²/τ_m) δ(t-t')δ_ik δ_jl`. **This is
   the mechanism that drives spontaneous bulk T1s** — without it (and
   without division), there is no stochastic driver for bulk edges to
   explore short lengths at all, regardless of how fluid `p0` nominally
   makes the tissue.

**T1 rule:** kinematic + energy-gated — "if an edge length goes below a
small threshold length `L_T1`, an intercalation occurs... if it results in
a lower energy" (their own stated rule, verbatim in spirit).

**Modeling the wound itself** (this is the part most relevant to the
"expansion then closure" question):
- Any cell fully or partially within radius `R_w` is removed. Actomyosin
  cortices in wound cells are disrupted (**contractility removed**, not
  boosted).
- The wound polygon keeps a residual **area-elastic modulus `K_w`**, which
  **decays to zero over 10 minutes**. While `K_w` is still nonzero, it
  resists further shrinkage; once tissue pre-stress is released by cutting
  the cells, the still-elastic hole plus the surrounding tissue's own
  tension produces a **rapid expansion of wound area**, which then
  contracts back toward the original size over ~10 minutes as `K_w` decays
  and the purse-string ring takes over. This exact area-vs-time pattern
  (rise to 125-175% then fall) is shown for real wing-disc wounds across
  multiple supplementary figures (Supp. Figs. 1, 7-9): every single wound
  in the dataset, regardless of genotype, rises before it falls.
- Purse-string line tension on the wound margin ramps from the mean line
  tension `Λ̄0` up to `Λ̄ps` (their Supp. Fig. 2d shows this ramp completing
  within roughly the first ~10-15 min).

**Default calibrated parameters (their Supp. Table 2):**

| Parameter | Symbol | Value |
|---|---|---|
| Normalized contractility | `Γ̄` | 0.04 |
| Normalized cell line tension | `Λ̄0` | 0.00 |
| Normalized purse-string tension | `Λ̄ps` | 0.08 |
| Preferred cell area | `A0` | 16 μm² |
| Friction coefficient | `μ` | 30 s |
| Line tension deviation | `σ_m` | 0.025 |
| Tension recovery time | `τ_m` | 150 s |
| Wound radius | `R_w` | 5.33 μm |
| Timestep | `Δt` | 3 s |
| Intercalation threshold length | `L_T1` | 0.2 μm |
| Mean % cells dividing/hour | — | 2.25 |

**Genotype-specific fits (their Supp. Table 3, used for WT/Rok/Mbs
comparison):**

| Parameter | Rok RNAi | WT | Mbs RNAi |
|---|---|---|---|
| `Λ̄0` | -0.141 | 0.012 | 0.119 |
| `Λ̄ps` | 0.047 | 0.150 | 0.096 |
| `μ` | 80 s | 80 s | 80 s |
| Division rate (%/hr) | 2.240 | 2.274 | 0.8723 |

### Key results

- **Row-resolved intercalation rate**: row 1 (wound margin) intercalates at
  a dramatically higher rate than row 2+ — statistically significant in
  every condition tested (WT: KS test D=1, p=0.0079; high purse-string
  simulation: D=1, p<0.0001; low purse-string simulation: D=0.9167,
  p<0.0001, with row2 vs row3 also significant at p=0.0337 in the
  low-purse-string case). Order of magnitude: row 1 rates run
  ~0.02-0.10 intercalations/μm²/hr/min depending on the exact metric and
  condition; row 2+ drops sharply, roughly 3-5x lower.
- **Enabling intercalation lowers cell energy at the wound margin** (Supp.
  Fig. 10): with intercalations enabled, time-averaged mechanical energy in
  the first row is measurably *lower* than with intercalations disabled
  (unpaired t-test, t=14.8, df=22, p<1e-12) — intercalation isn't just a
  side effect of stress, it's an energetically favorable relief valve.
- **Intercalation is necessary for closure in their own model**: Video
  legends 3 and 4 state this explicitly — "Vertex model simulation of wound
  healing with intercalations disabled: ... **the wound fails to close**
  during the simulation" vs. "with intercalations enabled ... **the wound
  is able to close**," same parameters otherwise.
- **Mbs RNAi (more solid/jammed) wounds fail to close** in real tissue too
  (Supp. Fig. 6a: "NA" closure time, i.e., no closure observed) — matching
  the simulation's solid-limit prediction. Rok RNAi and ubi-Ecad-GFP wounds
  close *faster* than WT (both p<0.0001/p=0.0001 vs WT).
- **Wound-area vs. junction-count trajectories differ by genotype**:
  Mbs RNAi wounds (Supp. Fig. 7) frequently show wound area fluctuating or
  plateauing without junction loss — consistent with a jammed margin that
  can't intercalate away stress. Rok RNAi and WT (Supp. Figs. 8-9) show
  area and junction count declining together, roughly in step.
- Closure time increases monotonically with contractility `Γ̄` and
  decreases with line-tension noise `σ_m` (Supp. Fig. 4a-d) — more
  fluctuation-driven fluidity closes wounds faster, more contractility
  alone (without fluidity) slows it down.

### Relevance to this codebase
This is the primary source for: the fluid/solid `p0` framework already used
(`SingleCellContractionParameters.m`), the row1≫row2+ intercalation
enhancement (approximated in this repo via `T1_TOL_wound` >
`T1_TOL_bulk`, since the bulk line-tension-fluctuation mechanism that
produces it mechanistically in their model is **not yet implemented here**
— flagged as an open gap), the energy-gate rule on T1 (`checkT1transitions.m`),
and the **missing wound-expansion-then-closure physics** (`K_w` decay term)
that is the most direct, still-unimplemented piece of this paper's model.

---

## 3. Kaurin, Bal & Arroyo 2022 — "Peeling dynamics of fluid membranes bridged by molecular bonds: moving or breaking"

**Citation:** Kaurin D, Bal PK, Arroyo M. *J. R. Soc. Interface* 19, 20220183
(2022). `rsif.2022.0183.pdf`

### What problem they solve
Not a tissue/vertex model at all — a continuum, reaction-diffusion-mechanics
model of **two adhered vesicles peeling apart**, where the adhesive patch
between them is made of mobile molecular bonds (cadherin-like) embedded in
a fluid membrane. The point: unlike solid-solid adhesion (fixed receptor
positions), bonds here can *move laterally* within the membrane in addition
to breaking and reforming — a genuinely different, richer failure physics.

### The model
State variables: bond concentration `c1(s,t)` and free-binder concentration
`c2(s,t)` on the patch, coupled to membrane mechanics (bending + capillary
tension) via a micromechanical model resolving the separation profile
`h(x,t)`.

- **Reaction-diffusion equations** for bonds/binders:
  `ċ1 = D1[c1' + c1(h'/x_γ²)']' + kon·c2² - koff·c1`
- **Force-dependent binding/unbinding rates (Bell's model):**
  - `kon(h) = k̄on · exp(-(h/x_γ)²)` — binding rate falls off as
    separation grows.
  - **`koff(h) = k̄off · exp(kh/f_β) = k̄off · exp(h/x_β)`** — this is
    Bell's law: unbinding rate rises *exponentially* with the force
    pulling the bond apart. `f_β` is the **force sensitivity** of the bond
    (a molecular property); `x_β = f_β/k` is the corresponding separation
    sensitivity.
  - Reference force scale for order-of-magnitude comparison:
    `f_γ = √(k·kB·T) ≈ 1 pN` for typical cadherin-like bond stiffness
    (`k ≈ 2.5×10⁻⁴ N/m`) at room temperature — **this is a real, dimensioned
    number in their vesicle system and does not transfer to a dimensionless
    vertex-model energy scale** (this is why this repo's `f_beta_T1`/`f_beta_T2`
    had to be calibrated empirically from this model's own force output
    rather than taken from this paper directly).

### Three distinct dynamical regimes, depending on the balance of
diffusion vs. reaction (Damköhler number) and on bond mobility:

1. **Diffusion-dominated** (long-lived, mobile bonds): the problem reduces
   to a classical Stefan (moving-boundary) problem. As force increases,
   bonds diffuse toward the shrinking patch and *concentrate*, which
   **self-stabilizes** the adhesion — the effective fracture energy
   increases as the patch shrinks. Similarity solution collapses onto a
   single universal curve (`X(τ) = 2λ√τ`).
2. **Reaction-dominated** (fast off-rate, low diffusivity — the regime
   appropriate for **cadherins**, since they're short-lived and partially
   immobilized by the actin cytoskeleton): decohesion proceeds via
   **traveling-front (FKPP-like) solutions** — a localized process zone at
   the patch edge moves at a constant, well-defined speed
   `v0 = √(D1·k̄off)`, with a small length scale `ℓ3 = √(D1/k̄off) ≈ 5 nm`
   controlling the front width. This is a **tear-out** mechanism, distinct
   from classical tear-out of immobile bonds because it still depends on
   small-scale diffusion near the moving front, not just marginal
   unbinding.
3. **Mixed reaction-diffusion regime**: multi-phasic, multi-timescale
   dynamics. Failure time follows a **power law**,
   **`t_fail ∝ (F - F_c)^(-2.2)`**, with a critical force `F_c` below which
   the patch is long-lived and above which it fails rapidly — `F_c` is
   interpreted as the mesoscopic strength of the adhesion patch, itself set
   by the microscopic bond force-sensitivity `f_β`. Force sensitivity of
   slip bonds (`f_β < f_γ`) sharply reduces this strength.

### Effect of crowding
With a finite maximum bond concentration `c_max` (rather than the dilute
limit), crowding **weakens** the patch and accelerates failure relative to
the dilute limit — the opposite of the diffusion-dominated
self-stabilization above. This is explicitly a *different* physical limit
(crowding vs. dilute-diffusion concentration) and the two should not be
conflated.

### Their own stated biological interpretation (important nuance)
They explicitly note that diffusion-dominated self-stabilization is the
right physical picture for **long-lived, freely mobile bonds**, but that
**cadherin junctions — anchored to the cytoskeleton, short-lived — are
better described by the reaction-dominated/tear-out limit**, which is the
one that does *not* self-stabilize and instead exhibits the accelerating
`(F-Fc)^-2.2` runaway failure.

### Relevance to this codebase
This is the direct source of the Bell's-law stochastic transition formula
now implemented in `checkT1transitions.m` and `checkWoundIntercalations.m`:

```
k_off = k_off0 · exp(F_drive / f_beta_eff)
P_flip = 1 - exp(-Δt · k_off)
```

and of the choice to make `f_beta_eff` **shrink as the edge shortens**
(`f_beta_eff = max(f_beta - slope·(tol - edge_len), 0.01)`) — a
phenomenological proxy for their reaction-dominated/tear-out regime
(the biologically appropriate one for cadherins, per their own framing),
**not** a literal implementation of their full reaction-diffusion bond
state `c1(s,t)` (which would require resolving a genuine PDE per edge —
out of scope here, and explicitly not what this repo claims to do).

---

## 4. How the three combine into this repo's strategy

| Paper concept | This repo's implementation | File(s) |
|---|---|---|
| Trepat: early crawling (OPTL) → later purse-string (IPTL) handoff | `crawl_force0` (active from t=0, decays) + `lambda_purse_string` (ramps in) | `getVertexForcesClassical.m`, `SingleCellContractionParameters.m` |
| Tetley: energy functional, `p0` jamming control | `E = ka(A-A0)² + (1/2rstiff)(P-p0)²`, `p0=4.0` for fluid regime | `getTissueEnergyClassical.m`, `SingleCellContractionParameters.m` |
| Tetley: T1 rule = length threshold + energy gate | `checkT1transitions.m` energy-gate block | `checkT1transitions.m` |
| Tetley: row1≫row2+ intercalation enhancement | `T1_TOL_wound` (0.15) vs `T1_TOL_bulk` (0.05) — a proxy, not the mechanism | `SingleCellContractionParameters.m` |
| Tetley: wound `K_w` decay → expansion-then-closure | **Not implemented** — flagged open gap | — |
| Tetley: bulk myosin line-tension fluctuations (fluidization source) | **Not implemented** — flagged open gap | — |
| Kaurin-Arroyo: Bell's-law stochastic gate | `k_off`/`P_flip` formulas, calibrated `f_beta_T1`/`f_beta_T2` | `checkT1transitions.m`, `checkWoundIntercalations.m` |
| Kaurin-Arroyo: reaction-dominated/tear-out (shrinking `f_beta_eff`) | `f_beta_length_slope_T1`/`_T2` | same |
| Kaurin-Arroyo: T1 (tangential) vs T2 (normal) force distinction | `total_tension` (T1) vs `pressure_area` (T2), never conflated | `getEdgeForces.m` |

### Open gaps, honestly, as of this summary
1. **Wound expansion phase** (Tetley's `K_w` decay) is the most direct,
   still-unimplemented piece — this repo's wound only ever contracts, it
   never gets the initial elastic-recoil expansion every real wound in
   Tetley's dataset shows.
2. **Bulk line-tension fluctuations** (Tetley's actual fluidization
   mechanism) are absent — `p0` alone shapes the energy landscape but
   provides no stochastic driver for bulk cells to explore short edges,
   so the row1≫row2+ falloff is currently approximated by a tolerance
   ratio rather than emerging mechanistically.
3. **`k_off0`** (the overall pace of the Bell's-law gate) has no
   paper-derivable value and has not been independently calibrated —
   `f_beta_T1`/`f_beta_T2` were calibrated from this model's own logged
   forces, but the base rate constant was not.
4. Trepat's crawling force is a pattern-level analogy (no substrate in
   this model), not a literal reproduction of traction magnitudes.
