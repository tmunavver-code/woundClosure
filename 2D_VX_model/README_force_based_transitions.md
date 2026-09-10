# Force-Based T1/T2 Transition Strategy — Implementation Spec

This document specifies how to replace the purely *kinematic* topological-transition
rule in the 2D vertex model (`2D_VX_model`) with a *force-based* one. It is written as
an implementation guide: each phase lists the goal, the files touched, the concrete
changes, and the acceptance check that tells you the phase is done.

The scientific goal: decide whether a cell edge undergoes a **T1 transition**
(neighbor exchange / intercalation) or a **T2 transition** (cell extrusion / removal)
based on the **traction on the edge**, not on an arbitrary length or area threshold.
A single force-dependent rule should reproduce both solid-like (rare intercalation,
purse-string closure) and fluid-like (intercalation-driven closure) behavior. That
unification is the novel contribution.

---

## 0. Background: what the current code does

- **T1 rule (`checkT1transitions.m`)**: flips an edge when `edge_length < T1_TOL`
  (a geometric threshold), with a separate, more permissive tolerance at the wound
  margin. No force enters the decision.
- **T2 / extrusion (`checkWoundIntercalations.m`)**: handles 3-cell extrusion on a
  geometric basis. `2D_VX_model` uses the corrected version; `2D_VX_model_reverted`
  has a known bug where extruded cells stay tethered to a master vertex — do not use it.
- **Force computation (`getEdgeForces.m`)**: already computes per-edge tangential
  (line-tension) and normal (area-pressure) forces, plus force-per-unit-length. This
  is the measurement layer we build the decision on.

**Physical mapping (fixed convention for this whole strategy):**
- T1 (edge collapse) is driven by the **tangential** line tension → gate on
  `force_tangent_per_length`.
- T2 (cell separation / extrusion) is a **normal peeling** event → gate on
  `force_normal_per_length`.

**Energy-function mapping** (our `(P - P0)^2` form ↔ the Farhadifar/Tetley
`Γ·P² + Λ·L` form): expanding ours gives contractility `Γ = 1/(2·r_stiff)` and
effective line tension `Λ = -2Γ·P0`. Therefore **`P0` is the target shape index `p0`**,
and the fluid–solid jamming line sits near `p0* ≈ 3.81`. This lets the force threshold
reference proximity to jamming (Phase 4).

---

## Phase 1 — Instrument and validate the force field (DO THIS FIRST)

**Goal:** confirm the force field is trustworthy *before* any transition rule depends
on it. Do not change the transition logic yet.

**Files:** `runEdgeForceAnalysis.m`, `checkT1transitions.m` (logging only),
`validateTrepatForces.m`.

**Changes:**
1. Log, at every candidate T1, the tuple
   `[tstep, v1, v2, length_at_flip, tangential_tension_per_length, normal_pressure]`
   into `T1ForceLog`. (Already wired in `checkT1transitions.m`.)
2. Produce two diagnostic figures:
   - **Force-vs-length scatter** at the moment of each T1. This answers the core
     question: *is the current length threshold already an implicit force threshold?*
     (If short edges always coincide with high tangential force, length is a hidden
     force proxy; if not, a force rule is genuinely needed.)
   - **Trepat spatial-force profile** (`validateTrepatForces.m`): radial and azimuthal
     edge-force components binned by distance from the wound center.

**Correctness requirements for `validateTrepatForces.m`** (these were the bugs that
made the first profile pure noise):
- **Minimum-length guard:** exclude edges shorter than `0.5 * T1_TOL_bulk` before
  computing `F/L`. Near-degenerate edges make force-per-length blow up and produce
  phantom spikes. A `< 1e-10` skip is NOT enough.
- **Time-average, don't single-snapshot:** pool edges across many stored snapshots
  before binning. A single snapshot in a 15×15 box has too few edges per bin.
- **Adaptive bin count:** target enough edges per bin; cap between ~6 and 15 bins.
  Do not hardcode 30 bins.
- **Blank under-sampled bins as `NaN`** (do not plot them as 0) so the line breaks
  instead of drawing fake zero-crossings.
- **Show per-bin edge counts** (a companion bar panel) and **SEM error bars** so noise
  is visually distinguishable from signal.

**Acceptance check:**
- Time-averaged Trepat profile shows a *coherent* azimuthal band behind the wound edge
  (not isolated spikes), with most bins above the minimum-count line. If bins stay
  under-sampled even after time-averaging, the domain is too small for a radial profile
  — in that case the T1 force-vs-length scatter becomes the primary validation instead.
- Sign convention verified: radial unit vector points *outward* from wound center;
  positive radial = outward.

---

## Phase 1b — Fix logging/architecture bugs (blockers for trustworthy data)

These do not change physics but corrupt the data if left unfixed.

1. **Vertex-ID sort convention:** the T1 logger looks up edges via
   `find([edgeForceData.v1]==min(v) & [edgeForceData.v2]==max(v))`. Confirm
   `getEdgeForces.m` stores `v1 = min(vertID)`, `v2 = max(vertID)` for every edge.
   If not, transitions are silently dropped from the log.
2. **Stale forces within a sweep:** `checkT1transitions` may perform multiple T1s in one
   call, but `edgeForceData` is computed once per timestep. After the first swap,
   forces are stale for later edges in the same sweep. Acceptable for the Phase-1
   diagnostic; must be fixed before force-*gating* (recompute forces after each swap,
   or process at most one swap per edge per step).
3. **Two intercalation systems:** `checkWoundIntercalations` and `checkT1transitions`
   can both fire on wound-edge edges in the same step. Confirm they cannot double-flip
   the same edge. Decide which system owns wound-edge T1s.
4. **Redundant `getEdgeForces` calls:** it is currently called ~3×/timestep. Compute
   once per step and reuse, so logged forces are internally consistent.

**Acceptance check:** number of logged T1s equals number of T1s actually performed
(no silent drops); no edge is flipped by both systems in one step.

---

## Phase 2 — Force-based transition criterion (the core change)

**Goal:** replace the hard geometric threshold with a stochastic, rate-dependent,
force-driven rule modeled on Kaurin–Arroyo (2022) bond-breaking kinetics.

**Files:** `checkT1transitions.m` (or a new `checkForceTransitions.m`),
`getVertexForcesClassical.m` / `getEdgeForces.m` (state variable), `param` setup.

**Changes:**
1. **Per-edge dynamical adhesion strength (Kaurin–Arroyo insight):** the junction's
   resistance to separation is a *dynamical variable*, not a material constant.
   Give each edge a state `S_ij` (analog of bond concentration `c1`), initialized from
   line tension, that strengthens as the edge shortens (crowding) and weakens as it is
   pulled apart. Store and evolve it per timestep.
2. **Bell's-law firing probability (replaces the hard cutoff):** instead of
   `if length < T1_TOL → flip`, compute a per-step transition probability that rises
   exponentially with the driving traction:

   ```
   P_flip = 1 - exp( -dt * k_off )
   k_off  = k_off0 * exp( F_drive / f_beta )
   ```

   - For **T1**: `F_drive = force_tangent_per_length` (edge-collapse driving force).
   - For **T2**: `F_drive = force_normal_per_length` (peeling / separation force).
   - `f_beta` is the force-sensitivity scale (Kaurin–Arroyo anchor: characteristic
     bond force `f_gamma = sqrt(k · kB·T) ≈ 1 pN` as an order-of-magnitude reference;
     in vertex-model units it is a fit parameter).
3. **Length/crowding modulation (Gap 5):** make the threshold scale `f_beta` *decrease*
   as the edge shortens, so the transition self-accelerates near collapse — matching
   the Kaurin–Arroyo runaway `t_fail ∝ (F − F_c)^-2.2`.
4. **Keep the energy gate (Tetley–Banerjee):** only perform the (stochastically chosen)
   flip if it also lowers system energy. Force *gates* the transition; energy still
   guards consistency. Do NOT remove the energy check.

**Parameters to add to `param`:** `k_off0`, `f_beta`, `f_beta_length_slope` (crowding),
and a master switch `param.useForceBasedTransitions` so the old kinematic rule can be
toggled for comparison.

**Validated against Kaurin–Arroyo (2022, J. R. Soc. Interface 19:20220183):** their model
resolves a genuine reaction–diffusion state, bond concentration `c1(s,t)`, which can
either *self-stabilize* (diffusion-dominated regime: bonds concentrate as the patch
shrinks, increasing effective fracture energy) or *runaway-fail* (reaction-dominated /
"tear-out" regime, appropriate for short-lived, cytoskeleton-anchored bonds like
cadherins — their §3.2). Our `f_beta_eff` shrinking with edge length is a phenomenological
proxy for the reaction-dominated/tear-out limit only — we do not resolve a true diffusing
`c1(s,t)` and therefore cannot reproduce their diffusion-dominated self-stabilizing regime.
This is the biologically appropriate limit for cadherin junctions (their own choice for the
cadherin case), so the simplification is defensible, but it should be stated explicitly
rather than presented as the general mechanism. The `t_fail ∝ (F−Fc)^-2.2` runaway cited
above is their fig. 5f fit in exactly this regime.

**Acceptance check:** with `useForceBasedTransitions = false` the code reproduces the
old behavior exactly; with it `true`, transition timing correlates with logged edge
traction (high-traction edges flip preferentially), and the run remains numerically
stable (no divide-by-zero on collapsed cells — keep the `< 3 vertices` guard clause).

---

## Phase 3 — Calibrate against the two limiting regimes

**Goal:** show the one rule reproduces both ends of the epithelial spectrum.

**Fluid limit (Tetley–Banerjee 2019):** high shape index `p0` / high line-tension
fluctuation. Expected: high intercalation rate concentrated in wound row 1, dropping
sharply in rows 2+, and the wound *closes*. Reproduce their row-resolved intercalation
falloff and their finding that enabling intercalations *lowers* wound-adjacent cell
energy.

**Solid limit (Mbs-RNAi analog):** low `p0` / high uniform tension. Expected:
intercalations rare even under stress, and the wound *fails to close* (matches the
Tetley Mbs-RNAi result and the `_reverted` "no contractility increase → no closure"
observation).

**Force-field target (Trepat 2014):** the simulated radial vs. tangential traction
pattern around the wound should qualitatively resemble the OPTL/IPTL structure
(outward tractions at the leading edge, inward/tangential band just behind). This is a
pattern check on the force field, not a magnitude match — there is no substrate here.

**Acceptance check:** both regimes reproduced by changing only `p0` (and tension
fluctuation), with `f_beta`, `k_off0` held fixed.

---

## Phase 4 — Unification (the publishable claim)

**Goal:** demonstrate that a single force-dependent T1/T2 criterion, with one or two fit
parameters, continuously interpolates between solid and fluid closure — no separate
mechanisms, no per-regime hand-tuned geometric thresholds.

**Change:** tie `f_beta` to proximity to jamming: near-solid tissue (`p0 < p0*`) requires
higher traction to flip; near-fluid (`p0 > p0* ≈ 3.81`) flips easily. Sweep `p0` across
the transition and show a smooth crossover in intercalation rate and closure success.

**Acceptance check:** a single parameter sweep in `p0` produces a monotonic
solid→fluid crossover in closure behavior under one unchanged rule.

---

## File-by-file summary of required edits

| File | Phase | Edit |
|------|-------|------|
| `validateTrepatForces.m` | 1 | Min-length guard, time-average, adaptive bins, NaN empty bins, per-bin counts, SEM error bars |
| `getEdgeForces.m` | 1b, 2 | Verify `v1=min,v2=max` sort; add per-edge strength state `S_ij` |
| `checkT1transitions.m` | 1, 1b, 2 | Keep force logging; add Bell's-law `P_flip`; keep energy gate; add `useForceBasedTransitions` switch |
| `checkWoundIntercalations.m` | 1b, 2 | Resolve double-fire with T1 system; add normal-force gate for T2/extrusion |
| `getVertexForcesClassical.m` | 2 | Keep `<3 vertices` guard; evolve `S_ij` |
| `runEdgeForceAnalysis.m` | 1, 1b | Call `getEdgeForces` once/step; collect more snapshots; call `validateTrepatForces(..., 'timeAverage', true)` |
| `param` setup | 2 | Add `k_off0`, `f_beta`, `f_beta_length_slope`, `useForceBasedTransitions` |

---

## Do-not-break invariants

- Never gate a transition on a force field that has not passed Phase 1 validation.
- Keep the degenerate-cell guard clause (`< 3 vertices` skip) in force calc and plotting.
- Keep the energy check as a gate; force decides *whether to attempt*, energy decides
  *whether it is allowed*.
- Preserve the old kinematic rule behind a toggle so every force result has a baseline
  to compare against.
- T1 gates on tangential force; T2 gates on normal force. Do not conflate the two.

---

## Phase 5 — Active force terms beyond T1/T2 (validated against Tetley 2019 and Trepat 2014)

The phases above only replace the *topological transition rule*. Two more mechanical
ingredients were checked against the wound-healing literature and one gap was closed:

**Purse-string tension — validated, matches Tetley SI exactly.** Tetley et al. (Nat. Phys.
2019, SI eq. for `E`) use the identical energy functional
`E = Σ K(A-A0)²/2 + Σ Γ P²/2 + Σ Λ_ij L_ij`, ramp the wound-margin line tension from `Λ0` to
`Λps` over ~10 min post-ablation (SI Fig. 2d), and — critically — state their T1 rule as
"if an edge length goes below `L_T1`, an intercalation occurs... *if it results in a lower
energy*" (SI p.6). That is precisely our energy-gate in `checkT1transitions.m`; the force-
based criterion in this doc extends their rule rather than replacing its logic. Our
`lambda_purse_string` ramp (`getVertexForcesClassical.m`) implements the same ramp
mechanism, and the line-tension force `F = Λ·(unit vector v1→v2)` is the correct negative
gradient of `Λ·L_ij`.

**Wound-directed crawling — was MISSING, now added.** Trepat et al. (Nat. Phys. 2014,
"Forces driving epithelial wound healing") show wound closure has two temporally distinct
mechanisms: an early outward-pointing traction layer (OPTL) from lamellipodial crawling of
row-1 cells into the wound, and a later inward/tangential layer (IPTL) from purse-string
tension transmitted to the substrate. Before this change, the model only had the
purse-string/contractility side (`lambda_purse_string`, `ka_wound_factor`,
`contractility_wound`) — there was no active, independent crawling force, i.e. the
*kinematic-to-dynamical* gap. `getVertexForcesClassical.m` now
applies `param.crawl_force0 * exp(-tstep/param.tau_decay_crawl)` per wound-margin cell,
directed from that cell's barycenter toward the wound centroid, active from `t=0` and
decaying on roughly the timescale that `lambda_purse_string` ramps up — reproducing the
early-crawl/late-purse-string handoff. Toggle: `param.enableWoundCrawling`.
Caveat (already correctly noted in `validateTrepatForces.m`): there is no substrate in this
model, so this is a pattern-level analogy (two temporally offset active mechanisms), not a
literal reproduction of traction-force magnitudes or the tangential-focal-adhesion
mechanism Trepat identifies for the IPTL.

**Still missing — myosin line-tension fluctuations (Tetley's fluidization source).** Tetley
et al. drive spontaneous bulk T1s (their "tissue fluidity") with an Ornstein–Uhlenbeck noise
term on every edge's line tension, `dΛ_ij/dt = -(Λ_ij - Λ0)/τ_m + ξ_ij(t)` (SI, "Line
tension fluctuations"), on top of cell division. Our model has neither: `p0` alone shapes
the energy landscape (correctly, per the `p0* ≈ 3.81` jamming mapping already in this doc)
but does not by itself generate bulk activity — without fluctuations or division, bulk
edges have no stochastic driver to explore short lengths, so intercalation may stay
artificially confined to the actively contracting wound margin regardless of `p0`. This
would matter most for reproducing Phase 3/4's row-resolved intercalation falloff and the
solid→fluid closure crossover. Not implemented here — flagged as the next candidate change,
pending a decision on scope (adds `τ_m`, `σ_m` to calibrate, and interacts with the
already-tuned wound-margin parameters).
