# Force-Based Wound Closure Strategy

This is a design document, not a change log — nothing here has been
implemented yet. It lays out how to move the *entire* wound-closure
pipeline from timer/geometry-driven to force-driven, using the three papers
summarized in `research_papers_summary.md` (Trepat 2014, Tetley 2019 SI,
Kaurin-Arroyo 2022).

---

## 0. Diagnosis: what's still kinematic today, piece by piece

One piece of this pipeline — the T1/T2 topological transitions — was
converted to genuine force-based (Bell's-law) gating earlier this session:
`checkT1transitions.m` and `checkWoundIntercalations.m` now compute real
per-edge forces via `getEdgeForces.m` and use them, with `f_beta_T1`/
`f_beta_T2` calibrated from this model's own logged force distributions.
That part is done.

Everything else in the closure pipeline is still driven by **fixed timers,
instant switches, or static geometric tolerances**, not by force:

| What | Current implementation | Where |
|---|---|---|
| Wound creation | Instant cell deletion, no elastic hole, no residual resistance | `createWound.m` |
| Purse-string tension ramp | `current_lambda = lambda_purse_string * min(1, tstep/t_ramp)` — a pure step-counter | `getVertexForcesClassical.m`, `getEdgeForces.m` |
| Crawl force decay | `current_crawl = crawl_force0 * exp(-tstep/tau_decay_crawl)` — a pure step-counter | `getVertexForcesClassical.m` |
| Wound-margin contractility boost | Instant `ka*2`, `rstiff/3` multiplier, on/off with no ramp, independent of the cell's actual stress state | `getVertexForcesClassical.m` |
| Row1≫row2+ intercalation enhancement | A fixed, larger length tolerance (`T1_TOL_wound` vs `T1_TOL_bulk`) — a geometric proxy, not a force effect | `SingleCellContractionParameters.m` |
| Bond/junction "strength" in the Bell's-law gate | `f_beta_eff` is a function of *instantaneous edge length only* — no memory of load history | `checkT1transitions.m`, `checkWoundIntercalations.m` |
| `k_off0` (base transition rate) | Fixed constant, never calibrated against anything | `SingleCellContractionParameters.m` |

None of these are wrong, exactly — they were reasonable placeholders — but
none of them respond to what the tissue is actually doing mechanically.
This strategy replaces each one with something driven by locally computed
force or an evolving state variable, grounded in a specific claim from one
of the three papers.

---

## 1. What "force-based" should mean here, per paper

- **Trepat**: the OPTL→IPTL handoff is a real physical handoff — crawling
  recedes as the ring's tension is actually transmitted, not on an
  independent clock. In this model (no substrate to transmit tension
  *to*), the closest honest analog is coupling the crawl force's decay to
  how much purse-string tension has *actually built up*, not to elapsed
  time.
- **Tetley**: wound expansion is a real mechanical consequence — released
  bulk pre-stress acting against a still-resistant, decaying wound-hole
  elasticity (`K_w`). It should fall out of a force balance, not be
  scripted. Bulk tissue fluidity also comes from a real, physical
  mechanism (myosin line-tension noise) — not merely from `p0` being
  numerically above the jamming point with nothing to actually shake bulk
  edges loose.
- **Kaurin-Arroyo**: their central, novel claim is that junction "strength"
  is a **dynamical state**, not a material constant — it evolves based on
  load history. The `f_beta_eff` proxy already in this codebase (shrinks
  with instantaneous edge length) is a snapshot, not a history — a
  genuinely force-based version needs a persistent per-edge state that
  integrates over time.

---

## 2. Phased strategy

### Phase F1 — Replace fixed timers with force/state-triggered ramps

Currently `current_lambda` and `current_crawl` are pure functions of
`tstep`. Replace both with feedback rules:

- **Purse-string buildup**: instead of ramping against elapsed steps, ramp
  against a real signal — recruit tension at a rate proportional to the
  margin's own tension deficit, e.g.
  `d(lambda_current)/dt = k_recruit * (lambda_purse_string - lambda_current)`,
  where the recruitment itself only switches on once the wound-hole
  elasticity `K_w` (Phase F2) has decayed past some fraction — i.e., the
  ring builds as the passive resistance fades, which is the actual
  causal order in Tetley's account, not an independent parallel timer.
- **Crawl force handoff**: decay the crawl force in proportion to how much
  of the purse-string's target tension has actually been recruited, e.g.
  `crawl_force(t) = crawl_force0 * max(0, 1 - lambda_current/lambda_purse_string)`.
  This makes the OPTL→IPTL handoff an actual relationship between the two
  mechanisms (Trepat's finding) instead of two independently scheduled
  exponentials that happen to overlap in the right ballpark.

### Phase F2 — Force-based wound creation and expansion (Tetley's `K_w`)

- Give the wound hole a residual area-elastic modulus `K_w`, initialized
  high at ablation, decaying via a relaxation ODE
  (`dK_w/dt = -K_w/tau_Kw`) rather than a hardcoded step-countdown.
  While `K_w>0`, the hole resists deformation; combined with the
  surrounding tissue's own released pre-stress, this should produce a
  genuine passive expansion phase, then let the purse-string (once
  recruited via Phase F1) take over and close it — reproducing the
  rise-then-fall wound-area trace every real wound in Tetley's dataset
  shows (their Supp. Figs. 1, 7-9), rather than scripting that shape in.
- **This is the most direct, still-missing piece from the earlier
  diagnosis in this repo's history** — everything else in this strategy
  is refinement; this one is a genuine missing mechanism.
- Acceptance check: instrument actual wound polygon area vs. time (not
  yet tracked anywhere in this codebase — only edge-count and per-cell
  area are currently logged) and confirm the rise-then-fall shape emerges
  from the force balance, not from a scripted schedule.

### Phase F3 — Bulk force-based fluidity (Tetley's line-tension noise)

- Add an Ornstein-Uhlenbeck noise term to the tangential line tension of
  **every edge, bulk included** (not just the wound margin):
  `dΛ_ij/dt = -(Λ_ij - Λ0)/tau_m + ξ_ij(t)`, feeding directly into the
  tension term already computed in `getEdgeForces.m`/
  `getVertexForcesClassical.m`.
- This is itself a stochastic *force* process, so it composes naturally
  with the Bell's-law gate already in place — `F_drive` now has a genuine
  fluctuating component from real active noise, not just from static
  cell-shape-derived tension.
- This directly replaces the current crude proxy for row1≫row2+
  intercalation enhancement (`T1_TOL_wound` vs `T1_TOL_bulk`, a fixed
  geometric tolerance difference) with the real mechanism: wound-margin
  cells intercalate more because they are under genuinely higher, more
  actively driven tension (purse-string + contractility boost), which the
  force-based T1 gate already responds to — once bulk noise exists too,
  the row-resolved falloff should emerge on its own rather than being
  hand-coded via a wider tolerance.

### Phase F4 — Dynamical per-edge junction strength (Kaurin-Arroyo, simplified)

- Add a persistent per-edge scalar state `S_ij` (not the full
  reaction-diffusion `c1(s,t)` — a lightweight ODE proxy, deliberately
  chosen in the reaction-dominated/tear-out regime Kaurin-Arroyo
  themselves identify as appropriate for cadherin-like junctions):
  `dS_ij/dt = -(S_ij - S0)/tau_S + alpha * max(0, -dL_ij/dt)` — strengthens
  as the edge actively shortens (crowding), decays back toward baseline
  otherwise.
- Replace the current instantaneous-length-only `f_beta_eff` with one
  driven by this evolving state: `f_beta_eff = f_beta_T1 * (S_ij / S0)`.
  This is what actually makes junction "strength" a dynamical variable
  with memory, rather than a snapshot function of current length — the
  specific novel claim of the Kaurin-Arroyo paper.

### Phase F5 — Force-triggered margin contractility (no more instant switch)

- Replace the current instant, unconditional `ka_wound_factor`/
  `contractility_wound` multiplier with a mechanotransduction-style
  feedback: the boost ramps in via the same kind of ODE as Phase F1's
  purse-string recruitment, scaled by the cell's own locally measured
  tension deficit (from `getEdgeForces.m`) rather than a flat switch
  applied identically to every margin cell regardless of its actual
  state. This also resolves the standing issue (flagged earlier this
  session) that the current instant boost contradicts Tetley's finding
  that margin cortices are disrupted first, not reinforced first.

### Phase F6 — Calibration loop (extends the approach already used)

- Every new parameter introduced above (`k_recruit`, `tau_Kw`, `tau_m`,
  `sigma_m`, `tau_S`, `alpha`) gets the same treatment `f_beta_T1`/
  `f_beta_T2` already got this session: no invented numbers, no
  paper-borrowed numbers (none transfer, per the discussion in
  `research_papers_summary.md`) — run the model, log the relevant
  quantity's real distribution, anchor the parameter to it (typically the
  median, for the same heavy-tail reasons already documented).
- Extend the acceptance checks already used (no NaN/collapse, energy
  trend sane, wound-edge count behavior) to include the wound-**area**
  trajectory, which is not currently instrumented anywhere and is the
  actual quantity "did it close" depends on.

---

## 3. Dependency order

- **F2** (wound `K_w` expansion) has no dependencies — do this first. It's
  also the one that most directly answers "does the model reproduce
  expansion-then-closure," independent of everything else.
- **F1** (timer→feedback ramps) is most meaningful once F2 exists, since
  the purse-string recruitment trigger is defined in terms of `K_w`
  decaying — natural to build together or immediately after.
- **F3** (bulk noise) is independent of F1/F2, can be done in parallel.
- **F4** (dynamical `S_ij`) builds on the "per-edge ODE state" pattern
  F3 establishes — natural to do after or alongside F3.
- **F5** (margin contractility feedback) reuses F1's feedback-loop
  pattern — do after F1.
- **F6** (calibration) is continuous, not a single gate — applies
  incrementally as each phase lands, the same way `f_beta_T1`/`f_beta_T2`
  calibration was done this session.

---

## 4. Recap: what's already force-based vs. what this adds

**Already done** (this session, prior to this document): T1 and T2
transition gating both use real per-edge forces from `getEdgeForces.m`
combined with a Bell's-law stochastic rule, calibrated against this
model's own force output, plus an energy gate on top.

**What this strategy adds**: extends "force-based" to the rest of the
pipeline — wound creation and expansion, the purse-string/crawl handoff,
margin contractility, bulk fluidity, and junction-strength dynamics — so
the whole system responds to locally computed forces and evolving state
rather than fixed timers, instant switches, or static geometric
tolerances standing in for force effects.
