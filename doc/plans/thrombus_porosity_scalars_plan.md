# Development plan: porous thrombus, passive scalars and reaction rules

Status: proposal, nothing implemented. Written against commit `3b81de5`
(branch `refactor`). File and line references are to that commit.

Revision 2: a node whose resistance reaches a threshold `alpha_wall` becomes a
bounce-back wall instead of a drag node (section 3.1.1, decisions D4, D8, D9).

Revision 3: the maintainer set the default `alpha_wall = 1e3` (D8); section 3.1
states what "explicit" and "implicit" drag mean in this plan.

Revision 4: transient binding of platelets before permanent binding (tethering
and translocation along the surface, then arrest), as work package 6 after the
others (section 3.6, 4.7, decisions D10 and D11).

This document is the brief for the agents that implement the work. It has four
parts:

1. what the code does today, including defects found while reading it (section 2);
2. the target model, with equations and units (section 3);
3. the architecture, with the decisions that the maintainer has to confirm (sections 1 and 4);
4. work packages with tasks, files, tests and acceptance criteria (sections 5 to 7).

Read `AGENTS.md` and `doc/user_guide/developer/` first. Everything there
applies (C++11 library code, collective MPI calls, bit-for-bit reproducibility,
output does not change state, units documented at every interface, docs and
`CHANGELOG.md` updated with behaviour).

Terms used below:

| Term | Meaning |
|---|---|
| solid fraction `phi_m` | Volume fraction of solid component `m` (bound platelets, fibrin, ...) in a lattice node, 0 to 1. Porosity is `1 - sum_m phi_m`. |
| permeability `K` | Darcy permeability of a node [m^2], computed from the solid fractions by a permeability law. |
| resistance `alpha` | Dimensionless drag rate per time step, `alpha = nu * dt / K`. Zero in free fluid. This is the number the fluid solver uses. |
| `alpha_wall` | Resistance at and above which a node is turned into a bounce-back wall. |
| geometry wall | A node that is a wall from the start (vessel wall, stent, ...). |
| clot wall | A node that became a bounce-back wall because its resistance reached `alpha_wall`. |
| pending node | A node with `alpha >= alpha_wall` that is not converted yet (it waits for the next wall update, or a cell vertex is still in it). |
| binding node | A node platelets can bind to: a wall binding site set by the case, or a node with enough deposited material (3.4). |
| anchor node | A non-wall, non-binding node with a binding node as face neighbour. It lies one lattice spacing from the binding node's centre, half a spacing outside the bounce-back surface. Vertices of binding cell types are attracted to it (3.6). |
| tethered | A vertex within the capture range of an active anchor node; a cell with at least a set number of such vertices. |
| scalar, species | A passive concentration field `c_k` transported by the fluid (thrombin, ADP, ...). |
| rule | A local function that turns concentrations, shear rate and solid fractions at a node into rates of change of solid fractions and concentrations. |

---

## 1. Decisions for the maintainer

Agents proceed with the default in each row unless the maintainer overrides it.
D1 has a measurement gate in work package 0; the default of D8 has a validation
gate in work package 1.

| # | Decision | Default (recommended) | Alternative | Why the default |
|---|---|---|---|---|
| D1 | Where the per-node resistance `alpha` lives and how the drag is applied | One extra external scalar in the fluid descriptor, read by a new collision dynamics that applies the drag implicitly. The new descriptor becomes the default `DESCRIPTOR` of all builds. | (b) same code, but only in a library variant `hemocell_porous`; (c) separate scalar field plus a force pass that adds the drag to the external force before `collideAndStream()`. | Velocity stays consistent everywhere it is computed (interpolation, output, statistics, scalar advection); no extra per-step pass; communicated and checkpointed with the lattice. Cost: +8 bytes per node (176 to 184) and incompatible lattice checkpoints. See 4.1. |
| D2 | Scalar transport solver | New finite-volume solver on `MultiScalarField3D<T>`, flux-limited advection. The existing CEPAC lattice is retired once the new solver covers `cases/CEPAC`. | Repair and generalise the CEPAC advection-diffusion lattice. | Physiological diffusivities are far below what a BGK advection-diffusion lattice can represent at HemoCell's `dx`, `dt` (see 3.2), and a D3Q19 lattice costs 22 doubles per node per species. |
| D3 | `SOLIDIFY_MECHANICS` compile-time variant | Fold it into the base library and make binding a run-time feature, if `sizeof(serializeValues_t)` does not change (the flag should fit into existing padding; verify). | Keep the variant and add a second test runner linked against it. | `hemocell_test` links only the base library, so binding cannot be tested today. Variants also multiply with D1(b). |
| D4 | Legacy binding mode (nodes become bounce-back walls) | No separate mode. A component with the permeability law `impermeable` puts its nodes above `alpha_wall`, which reproduces the legacy behaviour through the one conversion mechanism of 3.1.1. | Keep a second code path. | One mechanism to test; the legacy case is its limit. |
| D5 | How rules are defined | C++ interface (`hemo::ReactionRule`) registered by the case, plus a small set of built-in rules configured from XML. No run-time expression parser. | Expression parser in XML. | No new dependency; the same pattern as `CellMechanics`/`addCellType`. |
| D6 | Source location | New files in `helper/` (subsystems), `core/` (dynamics), `config/` (descriptor). | A new top-level directory; then extend the CMake globs and include directories and the table in `AGENTS.md`. | No build-system change. |
| D7 | Concentration unit | The solver treats concentrations as opaque numbers; documentation and examples use mol/m^3. Rate constants are given in SI and converted with `dt`. | Fixed unit. | Passive scalars are linear in `c`; only rules give the unit a meaning. |
| D8 | Value of `alpha_wall` | Configurable (`<wallResistance>`), default `1e3` (set by the maintainer). Validation V9 in WP1 has to show that the drag scheme passes up to that value; if it does not, stop and report. Fatal at set-up if a configured value exceeds the validated bound recorded in the code. | A lower value; `2` is the stability limit of an explicit drag force, which this plan does not use. | With the implicit drag of D1 no resistance is unstable, so the threshold is a modelling choice: where a node carries so little flow that a wall represents it better. At production parameters `alpha = 2` is `K = 5.5e-14 m^2`; a threshold there would make most clot material impermeable (table in 3.1.1). |
| D9 | Scalars in clot walls | Clot walls stay part of the scalar domain as diffusion-only nodes (no advection; per-species factor on `D`, default 1). Geometry walls stay impermeable. | Clot walls impermeable to scalars as well. | Agonists are produced on and inside the clot; an impermeable clot would trap or lose them. |
| D10 | Level of the transient-binding model (WP6) | Deterministic and without bond objects: attractive anchor nodes, plus an optional friction force on tethered vertices that stands for the mean effect of bonds forming and breaking. Both terms confirmed by the maintainer. | Discrete stochastic bonds (a spring per vertex to a fixed wall point, formed and broken with rates). | The default needs no per-bond state, no random numbers and no new communication. An attraction alone is conservative: it holds a platelet near the wall and can pin it, but it does not dissipate, so the friction term supplies the tunable slowing. Discrete bonds are the more faithful model and can replace the force law later behind the same interface (3.6). |
| D11 | Where the adhesion force is stored | Added to `sv.force_repulsion` in the same update as the repulsion, before the common cap. The arrest counter is a `uint16_t` in the padding of the particle record, if `sizeof(serializeValues_t)` does not change (verify). | A separate force member in the particle record (+24 bytes per vertex, new checkpoint layout). | No change of the particle layout or of checkpoints. Cost: adhesion and repulsion share one output field and one update interval. |

---

## 2. Current state

### 2.1 What exists

**Platelet "solidification".** Compiled only into `hemocell_solidify_mechanics`
(`SOLIDIFY_MECHANICS`, which adds `bool solidify` to the particle record,
`core/hemoCellParticle.h:63`). Enabled by `<parameters><enableSolidifyMechanics>`
(`config/config.cpp:224`) and per type by `HemoCell::enableSolidifyMechanics(name)`.
Each `solidifyTimescale` iterations, `HemoCell::iterate()` runs
(`core/hemoCell.cpp:419-425`):

1. `prepareSolidification()`: for each enabled type, `CellMechanics::solidifyMechanics()`.
   `PltSimpleModel::solidifyMechanics()` (`mechanics/pltSimpleModel.cpp:223-274`)
   takes every complete cell with a vertex flagged `sv.solidify`, finds its inner
   lattice nodes with `OctreeStructCell::findInnerNodes()`, turns each into a wall
   with `defineDynamics(..., new BounceBack)`, marks it as a binding site, and
   tags the cell's particles for removal.
2. `syncEnvelopes()`.
3. `solidifyCells()` (`core/hemoCellParticleField.cpp:931-984`): removes tagged
   particles, then for every binding site looks at the particles in the 3x3x3
   node neighbourhood and sets `sv.solidify` when the particle is within
   `distanceThreshold` and the local "Tresca" value exceeds `shearThreshold`
   (both from `<MaterialModel>` in the cell XML).

**Binding sites.** `helper/bindingField.*`: a `MultiScalarField3D<bool>` on the
domain lattice's block structure plus a `std::set<Dot3D> bindingSites` per
atomic block. Initial sites are wall nodes next to fluid inside a box
(`HemoCellFields::populateBindingSites()`, `core/hemoCellParticleField.cpp:840-867`).
Checkpointed as `bindingSites.{dat,plb}`; written as fluid output
`OUTPUT_BINDING_SITES` (`io/FluidHdf5IO.hh:312-335`).

**Passive scalar ("CEPAC field").** One `MultiBlockLattice3D<T, CEPAC_DESCRIPTOR>`
(advection-diffusion D3Q19), created when `<enableCEPACfield>` is set
(`core/hemoCellFields.cpp:118-138`), coupled with `LatticeToPassiveAdvDiff3D`
integrated into the fluid lattice, stepped right after the fluid
(`core/hemoCell.cpp:406-410`), checkpointed as `CEPAC.{dat,plb}`, written by
`writeCEPACField_HDF5()`. A block pointer is passed to `solidifyMechanics()` but
no model uses it. `cases/CEPAC` is the only user.

**Infrastructure that the new features touch.**

* `iterate()` order and the force reset at its end (`core/hemoCell.cpp:377-463`);
  the external force is zero between iterations unless the case sets a body force.
* Fluid dynamics: `GuoExternalForceBGKdynamics<T, DESCRIPTOR>` with
  `DESCRIPTOR = ForcedD3Q19Descriptor` (three external scalars: the force).
  The lattice is created by `HemoCell::initializeLattice()`,
  `initializeLatticeFromFlagMatrix()` or directly by about 30 drivers and tests.
  In Palabos the velocity is `u = j/rho + F/2` with `F` taken from the external
  field as an acceleration (`palabos/src/basicDynamics/externalForceDynamics.hh:68-80`).
* The fluid lattice is allocated with `<fluidEnvelope>` (2) but communicates one
  layer: the patched `changeEnvelopeWidth(1)` changes the overlaps only, and
  `getEnvelopeWidth()` keeps returning the allocated width. Auxiliary fields
  built from the lattice's management therefore get the same block-local
  indices as `atomicLattice` (the binding and interior-viscosity fields rely on
  this). A field with another envelope width needs `computeRelativeDisplacement()`.
* IBM kernels skip nodes whose dynamics `isBoundary()` and renormalise
  (`core/immersedBoundaryMethod.h:73-131`); `advanceParticles()` deletes
  vertices whose nearest node is a boundary node.
* Wall repulsion uses a node list built once (`populateBoundaryParticles()`).
* Fluid output: `WriteFluidField` (`io/FluidHdf5IO.hh`), one `case` per
  `OUTPUT_*` constant, datasets `(z, y, x, components)` in single precision.
* Checkpoints: `HemoCellFields::save/load` plus helper-specific files; the
  fingerprint covers the particle layout, cell types, `dx`, `dt`, not the
  lattice layout (`core/hemoCellFields.cpp:328-423`).
* Pre-inlet: two lattices; both multi-block objects exist on every rank, each
  rank owns blocks of one of them. Load balancing refuses to run with a
  pre-inlet or the CEPAC field (`helper/loadBalancer.cpp:40`).

### 2.2 Defects and limits found while reading

None of these were reproduced by a run; each follows from the code. They matter
because the new work either fixes them or must not copy them.

1. **The CEPAC field does not relax.** `param::tau_CEPAC` is initialised to 0 and
   never assigned (`mechanics/constantConversion.cpp`; the computation in
   `cases/CEPAC/CEPAC.cpp:55-59` is commented out), and it is passed as the
   `omega` argument of `AdvectionDiffusionBGKdynamics`
   (`core/hemoCellFields.cpp:127`). The value is also a relaxation time, not a
   rate.
2. **The CEPAC lattice uses the particle envelope width** (`envelopeSize`, 25
   in the shipped configurations) as its envelope (`core/hemoCellFields.cpp:124`).
3. **The shear criterion hard-codes `dt = 1e-7`:** `abs(tresca / 1e-7) > shearThreshold`
   (`core/hemoCellParticleField.cpp:973`). `tresca` is half the difference of the
   extreme eigenvalues of the strain-rate tensor, which is `gamma_dot / 2` in
   simple shear, so the threshold is not a shear rate in 1/s even at that `dt`.
4. **Binding flags can be lost at block boundaries.** `solidifyCells()` sets
   `sv.solidify` on any particle of the block, envelope copies included. The next
   `syncEnvelopes()` drops non-local particles and receives them again from
   their owners, so a flag set on a copy disappears. This matches the open item
   "Solidification communication" in `TODO.md`.
5. **Wall nodes created by binding are per block.** Each block converts the
   inner nodes that fall in its own bulk and envelope, and only if it holds the
   complete, flagged cell. Lattice dynamics in envelopes are not communicated.
6. **Cells that touch an aggregate are deleted.** New wall nodes are not added to
   the wall-repulsion list, so nothing keeps other cells away; a vertex that
   enters a converted node is removed and its cell is deleted as broken.
7. `solidifyTimescale` is not checked against the particle velocity update
   interval in `sanityCheck()`, although the step relies on fresh envelopes.
8. `distanceThreshold` is in lattice units in the cell XML (the convention is
   micrometres); `solidifyThreshold` is in the shipped `PLT.xml` files but is
   read nowhere.
9. The binding field's envelopes are never communicated; `refillBindingSites()`
   and the output read it by block.
10. The checkpoint fingerprint does not cover the lattice cell layout; a
    descriptor mismatch is only a warning (`BuildConfiguration`).

---

## 3. Target model

### 3.1 Porous medium (Darcy/Brinkman drag)

Momentum equation with a drag term in porous nodes:

    rho (du/dt + u . grad u) = -grad p + mu lap u + f_cells + f_body - (mu / K) u

In lattice units with Palabos' convention (external "force" is an acceleration,
`rho` close to 1), the drag is `F_D = -alpha u` with

    alpha = nu_lbm / K_lbm = nu * dt / K        (nu [m^2/s], dt [s], K [m^2])

`alpha` is dimensionless and spans many orders of magnitude: with
`nu = 1.1e-6`, `dt = 1e-7`, `K = 1e-14 m^2` gives `alpha = 11`, `K = 1e-17 m^2`
gives `alpha = 1.1e4`. An explicit drag force is unstable above `alpha = 2`, so
the drag is applied **implicitly**. With `G` the sum of all other forces on the
node (IBM spreading, body force) and `j` the momentum of the populations:

    u       = (j / rho + G / 2) / (1 + alpha / 2)
    F_total = G - alpha u

and the usual Guo collision is done with `u` and `F_total`. This satisfies
`u = j/rho + F_total/2`, is stable for every `alpha >= 0`, and tends to `u = 0`
for large `alpha`. For `alpha = 0` every operation reduces exactly to the
present dynamics (`x / 1.0` and `G - 0 * u` are exact), which gives a
bit-for-bit regression test.

This is the treatment of Guo and Zhao (Phys. Rev. E 66, 036304, 2002; check the
reference) without the porosity factor in the equilibrium and without the
Forchheimer term. In a uniform flow the scheme gives

    u_next = (u (1 - alpha/2) + G) / (1 + alpha/2)

so the velocity relaxes to `G / alpha` with the factor
`(1 - alpha/2) / (1 + alpha/2)` per step: its magnitude is below one for every
`alpha > 0`, and it is negative for `alpha > 2` (the population momentum
changes sign every step while `|u|` decays, which is what full-way bounce-back
does in the limit).

**Explicit and implicit drag.** The lattice Boltzmann step itself is explicit:
new populations follow from old ones, and no system of equations couples the
nodes. "Implicit" here concerns only how the drag term is evaluated within one
node. The velocity of a forced node is defined by `u = j/rho + F/2`. With a
drag, `F` depends on `u`, so this is an equation for `u`:

* *explicit*: evaluate the drag with the velocity without drag,
  `F_D = -alpha (j/rho + G/2)`. Per step the momentum is multiplied by
  `1 - alpha`: it oscillates for `alpha > 1`, diverges for `alpha > 2`, and the
  steady velocity is `G/alpha - G/2` instead of `G/alpha`;
* *implicit*: solve the equation. It is linear and involves this node only, so
  the solution is the closed form above, one division per node. No matrix, no
  iteration, no neighbour. The factor per step is
  `(1 - alpha/2) / (1 + alpha/2)` and the steady velocity is exact.

It is the same difference as between forward Euler and the trapezoidal rule for
`dy/dt = -a y`. What a lattice Boltzmann scheme cannot treat implicitly without
a global solve are terms that couple nodes (viscosity, pressure); a local term
can always be solved for in the node. The implicit formula can be evaluated
inside the collision (D1 default) or in a pass that writes the resulting force
to the external field (D1(c)); the two differ in where the result is stored,
not in stability.

The drag acts only below `alpha_wall`; at and above it the node becomes a
bounce-back wall (3.1.1). Pending nodes keep the drag until they are converted,
so the value used by the collision must stay finite (`inf * 0` is NaN in
`F_total`): cap it at an internal `alpha_max` (`1e9`).

**Solid components and permeability laws.** The state is a set of named solid
fractions `phi_m(x)`. `alpha` is derived, never stored as primary state:

    alpha = nu * dt * sum_m 1 / K_m(phi_m)      (terms with phi_m = 0 contribute 0; may be infinite)

    alpha <  alpha_wall : drag node with this alpha
    alpha >= alpha_wall : wall (pending until converted, with min(alpha, alpha_max) as drag)

Permeability laws to provide (per component, selected in XML):

| Law | `K(phi)` | Parameters |
|---|---|---|
| `linear` | `1/K = (phi / phi_ref) / K_ref` | `K_ref` [m^2], `phi_ref` |
| `kozenyCarman` | `K = d^2 (1 - phi)^3 / (180 phi^2)` | particle diameter `d` [m] |
| `davies` | `K = a^2 / (16 phi^1.5 (1 + 56 phi^3))` | fibre radius `a` [m] |
| `impermeable` | `K = 0` where `phi >= phi_ref`, no resistance below | `phi_ref`; gives a wall directly (legacy behaviour, D4) |
| custom | C++ functor registered by the case | |

Check the Kozeny-Carman and Davies forms and constants against the primary
sources before coding them. Literature permeabilities of fibrin gels and
platelet-rich clots span roughly 1e-13 to 1e-18 m^2 (e.g. Wufsus, Macera and
Neeves, Biophys. J. 2013); treat that as an order-of-magnitude guide and cite
verified values in the user documentation.

### 3.1.1 Switch to a bounce-back wall above `alpha_wall`

**What the threshold is.** The instability belongs to the *explicit* drag
force (drag computed from the velocity without drag): it diverges for
`alpha > 2`, and below that its steady velocity is off by a term of order
`alpha` (`G/alpha - G/2` instead of `G/alpha`). The implicit form of 3.1 costs
the same and has neither problem, so with the default of D1 no resistance is
unstable. The switch is adopted for these reasons instead:

* above some resistance the node carries a negligible flow and the Brinkman
  length `sqrt(nu_lbm / alpha)` is a small fraction of a node, so the drag
  resolves nothing that a wall does not;
* a wall is what the rest of HemoCell understands: IBM kernels skip it, wall
  repulsion acts on it, cells cannot be dragged into it;
* it bounds the regime in which the drag scheme has to be validated (V9).

`alpha_wall` is therefore a model parameter, default `1e3` (D8). It is far
above the explicit limit of 2, so the drag must be implicit.

What a threshold means at `nu = 1.1e-6 m^2/s`, `dt = 1e-7 s`, `dx = 0.5 um`
(`nu_lbm = 0.44`):

| `alpha` | `K = nu dt / alpha` | Brinkman length [nodes] | `1 / (1 + alpha/2)` |
|---|---|---|---|
| 2 | 5.5e-14 m^2 | 0.47 | 0.5 |
| 1e2 | 1.1e-15 m^2 | 0.066 | 0.02 |
| 1e3 | 1.1e-16 m^2 | 0.021 | 0.002 |
| 1e4 | 1.1e-17 m^2 | 0.0066 | 0.0002 |

Everything less permeable than the `K` of the chosen row is impermeable to
flow. With the range quoted above, a threshold of 2 would turn bound platelets
and much of the fibrin range into walls at once, which reduces feature 1 to the
present behaviour; `1e3` keeps Darcy flow through fibrin-like material and
leaves only dense platelet material as walls. Scalars still diffuse through
clot walls (D9).

**Rules of the switch.**

1. *State.* A stored flag field `clotWall` (one flag per node, envelope as the
   other porous fields). It changes only in the wall update below. It is
   checkpointed and not re-derived at restart, because the guard of rule 3
   makes it depend on where cells were.
2. *Wall update.* Runs every `wallUpdateEvery` iterations, a multiple of the
   particle velocity update interval, right after `syncEnvelopes()` (so every
   block holds complete cells), and after the binding step when both are due.
   Collective on every call, whether or not a node changes.
3. *Decision, by the block that owns the node (bulk only).* A node is converted
   if it is not a wall, `alpha >= alpha_wall`, and (with `protectCells`, default
   on) no vertex has it as nearest node. The particle grid gives the vertices
   per node. Particles of a cell that is bound in the same step are removed
   before the decision. Nodes that fail the guard stay pending and are counted
   (reduced over ranks, logged at output steps).
4. *Application.* Exchange the envelopes of `clotWall`; then every block gives
   every newly flagged node in its bulk **and envelope** bounce-back dynamics
   (`BounceBack<T, DESCRIPTOR>(1.)`, one object per node as `defineDynamics()`
   does) and sets the stored resistance of the node to 0. Envelope nodes must
   carry the same dynamics as their owners: Palabos collides the envelope layer
   locally before streaming it into the bulk (see the comment in
   `HemoCellParticleField::spreadParticleForce()`). The same holds for `alpha`.
5. *Populations* of a converted node are left as they are; bounce-back returns
   them in the next step. The node's fluid mass leaves the active fluid; sum
   and log the converted volume.
6. *One way.* In this version a clot wall stays a wall even if a rule lowers
   `phi` afterwards. Reverting needs new populations for the node and is an
   open question (section 8).
7. *Consumers that cache `isBoundary()`* are updated in the same call, on bulk
   and envelope: the wall-repulsion node list (a converted node next to a
   non-wall node is added), the binding node list, the anchor node list (3.6),
   the scalar mask (4.3) and the per-block fluid-node count. A unit test compares each cache with a
   rebuild from scratch.

The effective wall position of a clot wall is that of bounce-back (half-way to
the neighbouring fluid node), while a drag node just below the threshold blocks
the flow roughly at the node. Validation V4 measures the resulting jump in flow
rate when a layer converts.

### 3.2 Passive scalars

For each species `k`:

    dc_k/dt + u . grad c_k = div(D_k grad c_k) + S_k

`u` is the fluid velocity of the lattice (in porous nodes this is the Darcy
velocity). Geometry walls are impermeable and not part of the scalar domain.
Clot walls (3.1.1) remain part of it as diffusion-only nodes (D9): they keep
their concentration when they convert, rules keep acting in them, and species
diffuse through them with `D_k` times a per-species factor (default 1).

**Why not the advection-diffusion lattice.** Protein and small-molecule
diffusivities are of order 1e-11 to 1e-10 m^2/s. With `dx = 5e-7`, `dt = 1e-7`,
`D = 5e-11` gives `D_lbm = D dt / dx^2 = 2e-5`, a relaxation time of
`0.5 + 6e-5`, and a grid Peclet number of 500 at `u_lbm = 0.01`. A BGK
advection-diffusion lattice is not usable there; keeping it stable
(`tau >= 0.51`) means `D >= 8e-9 m^2/s`, two orders of magnitude too diffusive.
First-order upwind differencing is no better: its numerical diffusivity is
about `u/2` in lattice units, hundreds of times `D_lbm`.

**Scheme.** Explicit finite volume on the lattice nodes (node = cell centre):

* advection with a second-order flux-limited (TVD) face value, limiter
  selectable (`vanLeer` default, `superbee`, `minmod`, `upwind` for reference);
  face velocity is the mean of the two node velocities and zero on faces that
  touch a wall node;
* diffusion with central differences and a face diffusivity (constant per
  species now; the interface takes a per-face value so a porosity-dependent
  `D_eff` can be added);
* the update preserves a uniform field exactly: use
  `c_new = c - dts * (div(u c)_limited - c * div(u)_faces) + ...` with the same
  face velocities in both terms. The IBM velocity field is not discretely
  divergence-free, and the conservative form alone would create concentration
  peaks at moving membranes. Total scalar mass is then conserved only up to the
  discrete divergence; it is monitored (test S5) and reported in the log;
* stencil width 2, so scalar fields have an envelope of 2 and one envelope
  exchange per update.

**Time step.** The scalar step is `dts = N_s * dt`, run every `N_s` fluid steps
(`N_s >= 1`, from XML). Stability checks at set-up and, for velocity, at each
update on every block (reduced over ranks, fatal on violation with the
measured value):

    (|u_x| + |u_y| + |u_z|)_max * N_s <= CFL_max (default 0.3)
    D_lbm * N_s <= 1/6

**Boundary conditions.** Per domain face: `periodic` (must agree with the
lattice periodicity), `dirichlet <value>`, `zeroGradient` (outflow), `noFlux`.
Faces of geometry walls are always `noFlux`. Faces of clot walls carry no
advective flux (face velocity zero) and a diffusive flux as described above.

**Sources.** `S_k` is the sum of (a) rule outputs (3.3), (b) wall-flux patches:
a flux `J` [unit m^-2 s^-1] on the fluid faces of wall nodes inside a box or
flagged as binding sites (injury site), (c) volumetric box sources, for tests.

**Resolved, not modelled.** Do not add an "RBC-enhanced" effective diffusivity:
the velocity fluctuations caused by the resolved cells advect the scalar. The
scalar passes through membranes (as the CEPAC field does); state this as a
model limitation.

### 3.3 Reaction rules

Every `N_r` fluid steps, for every fluid node and every clot-wall node of the
domain (geometry walls are skipped):

    inputs : c_k (all species), phi_m (all components), shear rate gamma_dot [1/s], dtr = N_r * dt [s]
    outputs: dc_k/dt [unit/s], dphi_m/dt [1/s], summed over all rules

followed by an explicit Euler update, clamping (`0 <= phi_m`, `sum_m phi_m <= phi_max <= 1`,
`c_k >= 0` by limiting sinks to `c_k / dtr`), an envelope exchange of the changed
fields, and a refresh of `alpha`. The number of clamped nodes per update is
summed over ranks and logged at output steps, never per node.

Built-in rules (all parameters in SI, from XML):

| Rule | Rate | Use |
|---|---|---|
| `deposition` | `dphi_m/dt = k * g(c_s) * h(gamma_dot) * (phi_max - phi_m)`; optional consumption `dc_s/dt = -y * dphi_m/dt` | fibrin formation driven by thrombin |
| `release` | `dc_s/dt = k * phi_m * h(gamma_dot)` | agonist or thrombin generated on bound platelets |
| `decay` | `dc_s/dt = -k * c_s` | inhibition, degradation |
| `conversion` | `dc_a/dt = -k c_a`, `dc_b/dt = +y k c_a`, optionally gated by `g(c_s)` | simple cascades |

`g(c)` is `step` (`c >= c_thr`), `hill` (`c^n / (c^n + c_half^n)`) or `one`;
`h(gamma_dot)` is a window `gamma_min <= gamma_dot <= gamma_max` (defaults 0, infinity).

Cases add their own rules by deriving from `hemo::ReactionRule` (section 4.4).

**Shear rate.** One library function returns the scalar shear rate
`gamma_dot = sqrt(2 S:S)` in 1/s for a fluid node, `S` being the strain-rate
tensor. For simple shear this equals the velocity gradient. It replaces the
hard-coded expression of defect 3 and is used by rules, by binding and by
output. With Guo forcing the non-equilibrium stress contains a force term
(`-(uF + Fu)/2`); include the correction or compute `S` from velocity
differences, and validate both in Couette and Poiseuille flow (test R3). Inside
dense porous nodes the value is not meaningful; rules that gate on shear should
be read as acting at the clot surface. In clot-wall nodes the shear rate passed
to rules is 0.

A rule may push `alpha` past `alpha_wall`. The node then becomes pending and
keeps the drag until the next wall update (3.1.1); the reaction step itself
never changes lattice dynamics.

**Time scales.** Clot growth takes seconds to minutes; a cell-resolved run
covers milliseconds to seconds. Rate constants are plain user inputs; if a case
accelerates biochemistry it does so by scaling them, and says so. The library
does not hide an acceleration factor.

### 3.4 Binding with porous deposition

Feature 1 of the request: a platelet that binds is removed and its volume
becomes solid fraction of the component `platelet`. Its nodes turn into walls
only where the resulting resistance reaches `alpha_wall`.

Per binding step (interval `N_b`, a multiple of the particle velocity update
interval), all steps collective and in a fixed order:

1. **Flag (owners only).** Each block examines only particles in its bulk
   (`localDomain`) of binding-enabled types. A particle is flagged when a
   binding node lies within `bindingDistance` [um] and the shear rate at the
   particle's node is inside the configured window (binding mode `direct`; in
   the mode `tethered` of 3.6 the arrest criterion sets the flag instead, and
   the rest of this step is unchanged). Binding nodes are read from
   bulk and envelope of the binding-flag field, whose envelope is up to date.
2. **`syncEnvelopes()`** so that every copy carries the flag.
3. **Deposit.** Each block, for each complete cell with a flagged vertex:
   compute the inner nodes; for those in the block's **bulk** add
   `phi_bound` (`boundSolidFraction` in the cell XML) to the `platelet`
   component and clamp. Then remove all particles of the cell from the block. If a
   cell has no inner node on any block (thin, unfavourable orientation), the
   node nearest to its centroid receives the deposit; decide this from
   quantities that every block computes identically.
4. **Exchange envelopes** of the changed fields; refresh `alpha`, the binding
   flags and node lists (rule below) and the repulsion node lists on bulk and
   envelope.
5. **Wall update** (3.1.1) if due: deposited nodes whose resistance reached
   `alpha_wall` become clot walls. With the `impermeable` law this happens for
   every deposited node in the same step, which is the legacy behaviour (D4),
   now decided per node by its owner and applied on envelopes as well.

Every node is written by exactly one block, the decision travels with the
particles, and nothing depends on message arrival order. A bound-cell counter
and the deposited solid volume are accumulated (reduced over ranks), written to
the log and to a CSV file, and stored in checkpoints.

**Other cells near the aggregate.** Drag nodes are not walls: vertices are not
deleted there and IBM kernels include them. Clot walls are walls: a vertex
whose nearest node is one is deleted and its cell with it, as at any wall
(defect 6). Two things keep free cells from being lost:

* repulsion starts before the wall: nodes with `sum phi >= phi_repulsion` next
  to a node below that value are in the wall-repulsion node list when it is
  enabled (and leave it if a rule lowers `phi`); clot walls are in the list as
  walls. Choose `phi_repulsion` so that its resistance is below `alpha_wall`;
* the `protectCells` guard of 3.1.1 delays the conversion of a node while a
  vertex is in it.

Without repulsion, cells that enter dense drag nodes are slowed to the local
Darcy velocity, and pending nodes that they occupy do not convert; document
this and recommend enabling wall repulsion with the porous medium.

**Binding to deposited material.** A node is a binding node when its solid
fraction, summed over the components flagged `bindsPlatelets`, reaches
`phi_binding`, whether it is a drag node or a clot wall. Wall binding sites set
by the case (`populateBindingSites()`) stay binding nodes. The flag field is
recomputed from these two sources at every refresh.

### 3.5 Order within `iterate()`

    1  repulsion (wall list now includes dense drag nodes and clot walls) and adhesion (3.6, WP6)
    2  spread particle forces
    3  collideAndStream               (drag is inside the collision; clot walls bounce back)
    4  scalar transport, if due       (velocity field, advect, diffuse, sources, exchange)
    5  reaction rules, if due         (shear rate, rules, update c and phi, exchange, refresh alpha)
    6  interpolate velocity + syncEnvelopes, if due
    7  binding, if due                (flag owners, sync, deposit, remove, exchange, refresh)
    8  wall update, if due            (decide in bulk, exchange, bounce-back on bulk + envelope, caches)
    9  advance particles, constitutive model, interior viscosity, clean-up
    10 reset lattice force, iter++

Steps 4, 5, 7 and 8 are new or rewritten. `N_b` and `wallUpdateEvery` must be
multiples of the particle velocity update interval (check in `sanityCheck()`);
`N_s` and `N_r` are free, and `N_r` should be a multiple of `N_s` (warn
otherwise). Lattice dynamics change in step 8 only, so steps 4 and 5 always see
the walls that the last collision used.

### 3.6 Transient binding before permanent binding (WP6)

Goal: a platelet near a binding surface is first held by weak, transient bonds.
They slow it down while it can still move and roll along the surface. When it
has been slow enough for long enough, binding becomes permanent and the deposit
of 3.4 happens. This replaces the instant criterion of 3.4 step 1 for cell
types that select the binding mode `tethered`.

**Anchor nodes.** The attraction acts towards anchor nodes, not towards the
wall. An anchor node is a non-wall node that is not itself a binding node and
has a binding node among its six face neighbours. Its centre is one lattice
spacing from the binding node's centre and half a spacing outside the
bounce-back surface, so the attraction never pulls a vertex further than that
towards the wall. (A vertex whose nearest node is a wall node is deleted, so
an attraction towards the wall node itself would destroy the platelet.) The
anchor list is derived from the binding flags on bulk and envelope of each
block, in a fixed node order, and is rebuilt with the other caches whenever the
binding flags or walls change (3.1.1 rule 7). It therefore follows the clot:
platelets tether to deposited material as they do to the injury site.

An anchor is *active* when the shear rate at its node lies in the tether window
`[tetherShearMin, tetherShearMax]` (defaults 0 and infinity), evaluated at each
adhesion update with `hemo::shearRate()`. This lets a case restrict tethering
to high shear.

**Attraction.** For every vertex `x` of a tethering cell type and every active
anchor `a` with distance `d = |x_a - x| < r_t`:

    F_att = k_t * 4 (d / r_t) (1 - d / r_t) * (x_a - x) / d

`k_t` is the largest force [pN], reached at `d = r_t / 2`; the force is zero at
the anchor and at the capture range `r_t` [um], so it has no jump. The
contributions of all anchors in range are summed in list order. The neighbour
search is the one of wall repulsion (the 3x3x3 nodes around a list node), so
`r_t` is limited to one lattice spacing; larger values are a fatal error.

Wall repulsion keeps acting. With the usual cutoff above one lattice spacing it
still pushes outwards at the anchor, so the vertex rests where the two balance,
somewhat outside the anchor layer. Test T2 checks this position against the
force balance.

**Friction.** An attraction is conservative. It brings the platelet into the
slow fluid near the wall, and because the anchors form a grid, the potential is
corrugated along the surface: below some shear the platelet is pinned, above it
it slides. It cannot represent bonds that stretch and break, which is what
dissipates energy. The optional second term does that in the mean:

    F_fric = -gamma_t * w * v,      w = max over active anchors in range of (1 - d / r_t)

with `v` the vertex velocity (`sv.v`; binding surfaces are at rest) and
`gamma_t` in pN s/um. A population of `N` bonds per vertex with stiffness
`kappa` and lifetime `1 / k_off`, dragged at speed `v`, is stretched by
`v / k_off` on average and resists with `N kappa v / k_off`, so `gamma_t`
stands for `N kappa / k_off`. Check this reading against the literature on
bond-mediated friction before it goes into the user documentation. Friction
acts only on the vertices near the surface, so the shear flow on the rest of
the platelet turns it over those vertices: the platelet rolls or flips instead
of only sliding.

The friction force uses the velocity of the last update and is then held for
`N_a` fluid steps (the adhesion update interval), so it has a stability bound,
checked at set-up (fatal):

    gamma_lbm * N_a * n_v <= 0.5,     gamma_lbm = gamma_t * dt / dm

where `n_v` is the largest number of tethering vertices per lattice node,
estimated from the mean edge length of the cell mesh.

**Storage and cap.** Both terms are added to `sv.force_repulsion` in the same
update as the repulsion and before `capRepulsionForce()`, so `FORCE_LIMIT`
bounds repulsion and adhesion together (D11). Adhesion therefore uses the
repulsion update interval, which is already required to be a multiple of the
particle velocity update interval, and must be enabled together with wall
repulsion. Each block computes the force for all its particles, envelope
copies included, as repulsion does.

**Arrest.** Evaluated in the binding step (interval `N_b`), right after
`syncEnvelopes()`, by every block for every complete cell of a tethering type:

    tethered  : at least n_min vertices within r_t of an active anchor
    slow      : | mean of sv.v over all vertices | * dx / dt < arrestSpeed
    counter   : counter + 1 if tethered and slow, else 0     (saturating)
    permanent : counter * N_b * dt >= arrestTime  ->  set sv.solidify

The sums run over the vertices in vertex-ID order, so every block that holds
the complete cell obtains the same numbers bit for bit. Each block writes the
new counter and the flag to the vertices it owns; after the `syncEnvelopes()`
of 3.4 step 2 all copies agree, and steps 3 to 5 of 3.4 deposit the platelet.
The counter is a `uint16_t` per vertex in the particle record (D11); every
vertex of a cell carries the same value.

**What this model does and does not represent.** It gives capture near a
binding surface, slowing that can be tuned by `gamma_t`, motion along the
surface, and arrest followed by the deposit. It does not give stop-and-go
motion from individual bonds, a memory of where a bond formed, or
force-dependent rupture (catch or slip bonds). Lubrication between platelet and
wall is not resolved at this grid spacing, so the speed along the wall depends
on the stand-off set by `k_t` and the repulsion. Real translocation lasts
seconds; a case that covers milliseconds has to scale `arrestSpeed` and
`arrestTime` and must say so.

**Later replacement by discrete bonds.** The force law sits behind one function
per block (4.7). A model with one spring per vertex to a fixed point on the
surface, formed and broken with rates, would replace that function and needs
per-vertex state (anchor point, 24 bytes) and a random number generator keyed
by cell ID, vertex ID and iteration so that results do not depend on the rank
count. It is not part of WP6.

---

## 4. Architecture

### 4.1 Fluid: descriptor and dynamics (D1)

* `config/hemoDescriptors.h`: `hemo::descriptors::ForcedPorousD3Q19Descriptor<T>`,
  i.e. `plb::descriptors::D3Q19DescriptorBase<T>` plus an external-field struct
  with `numScalars = 4`, `forceBeginsAt = 0`, `sizeOfForce = 3`,
  `resistanceBeginsAt = 3`, `sizeOfResistance = 1` (follow
  `Force3dDescriptor`/`ForcedRhoBarJD3Q19Descriptor` in
  `palabos/src/latticeBoltzmann/` for the required members, `numSpecies` and
  the static `name`). It keeps the D3Q19 base, so Palabos' D3Q19 template
  specialisations still apply. `constant_defaults.h` selects it as `DESCRIPTOR`.
* `core/porousDynamics.h/.hh`: `hemo::GuoDarcyBGKdynamics<T, Descriptor>`,
  derived from `plb::ExternalForceDynamics`. It overrides `computeVelocity`,
  `computeVelocityExternal`, `collide`, `collideExternal`, `clone`, `getId`,
  and the serialisation constructor, following `GuoExternalForceBGKdynamics`
  (`palabos/src/basicDynamics/externalForceDynamics.hh:333-350`), and registers
  an id with `meta::registerGeneralDynamics` so that checkpoints can restore it.
  The Guo force term needs `F_total`, not the stored force: call the
  implementation template with a local 3-array or write the loop out.
* The lattice is created with this dynamics as **background dynamics** when the
  porous medium is enabled in the configuration: `initializeLattice()` and
  `initializeLatticeFromFlagMatrix()` choose it; drivers that create the lattice
  themselves use a small factory `hemo::newFluidDynamics(cfg)`. The
  `PorousMedium` constructor checks the background dynamics and stops with
  `fatal()` if it is the plain Guo dynamics. Do not assign the dynamics per node
  with `defineDynamics()` over a domain (that clones one object per node).
* Boundary-condition composites wrap a clone of the background dynamics; verify
  in WP0 that inlets and outlets behave as before. Interior viscosity clones it
  too and would inherit the drag, but the combination is refused for now (4.2).
* `setExternalVector(..., forceBeginsAt, ...)` resets three scalars, so `alpha`
  survives the force reset. `saveExternalForce()` in `io/FluidHdf5IO.cpp` copies
  all external scalars; that stays correct.
* Add the lattice layout (descriptor name, `numScalars`, `sizeof(T)`) to the
  checkpoint fingerprint and make a mismatch fatal with a clear message.

Rejected for D1: a pass that adds the drag force (implicit formula) to the
external field before `collideAndStream()`. The force field is reset at the end
of `iterate()`, so
velocities computed between iterations (output, statistics, any case code)
would miss the `F/2` term, which is not small in porous nodes (`u*` is
`(1 + alpha/2)` times the true velocity). Rejected as primary: per-node
dynamics objects carrying `alpha` (needs `dataStructure` communication for
envelopes and re-creation after restart; the interior-viscosity code shows how
fragile this is). It remains the fallback if a descriptor change is vetoed.

The wall switch of 3.1.1 does not replace the implicit drag. Pending nodes sit
above `alpha_wall` until the next wall update (longer with `protectCells`), and
the drag has to be exact and stable for them. With `alpha_wall = 1e3` (D8) the
explicit form of 3.1 is excluded: it would need `alpha_wall <= 1`, which makes
everything below `K = nu dt` (1.1e-13 m^2 at production parameters) a wall.

### 4.2 `hemo::PorousMedium` (`helper/porousMedium.*`)

Owned by `HemoCellFields` (created on demand, deleted in its destructor, like
`bindingHelper`). Built on the **domain** lattice's block structure with the
lattice's allocated envelope width, periodicity copied from the lattice.

    class PorousMedium {
    public:
      int addComponent(const std::string & name, PermeabilityLaw * law, bool bindsPlatelets);
      int component(const std::string & name) const;                 // fatal if unknown
      plb::MultiScalarField3D<T> & solidFraction(int component);
      plb::MultiScalarField3D<bool> & bindingNodes();                // replaces bindingFieldHelper's field
      void exchangeEnvelopes();          // collective; after any change of phi or binding flags
      void refreshResistance();          // collective; phi -> alpha on bulk + envelope, binding/repulsion lists
      plb::MultiScalarField3D<int> & clotWall();                     // 1 where a node was converted (3.1.1)
      void updateWalls();                // collective; convert pending nodes, call after syncEnvelopes()
      plint pendingNodes();              // collective reduction; nodes at or above alpha_wall not yet converted
      void checkpoint();  static void restore(HemoCellFields &);     // collective
      T totalSolidVolume(int component); // collective reduction, [m^3]
    };

Block-level code reaches the block's fields through pointers set on
`HemoCellParticleField` (as `bindingField` today) or through
`getComponent(atomicBlockId)`; never through `MultiScalarField3D::get()`, which
is collective.

`helper/bindingField.*` is absorbed: wall binding sites (`populateBindingSites`)
and deposited nodes share the one flag field and the per-block node list.

`refreshResistance()` never changes lattice dynamics; `updateWalls()` is the
only place that does (3.1.1, rules 2 to 7). Both write bulk and envelope of the
fluid block from fields whose envelopes were just exchanged, so no
`dataStructure` communication of the lattice is needed. `restore()` loads
`clotWall` and checks every local node: a node is a wall in the lattice exactly
if it is a geometry wall or flagged in `clotWall`, otherwise `fatal()`. The
geometry-wall mask is taken once, before the first conversion, and is
reproducible from the case set-up at restart.

Interior viscosity assigns dynamics to nodes inside cells
(`findInternalParticleGridPoints()`, `internalGridPointsMembrane()`) and resets
them to the background dynamics. Check that it never touches a wall node; until
a test shows it, the combination of interior viscosity and the porous medium is
a fatal error at set-up.

### 4.3 `hemo::ScalarTransport` (`helper/scalarTransport.*`)

Owned by `HemoCellFields`. Per species: name, `D` [m^2/s], initial value, six
face conditions, sources, one `MultiScalarField3D<T>` for `c` and one scratch
field shared by all species. One `MultiTensorField3D<T,3>` holds the velocity,
filled once per scalar step on bulk and first envelope layer from
`computeVelocity()` of the fluid block (the envelope layer is valid after
`collideAndStream()`; `interpolateFluidVelocity()` relies on the same). A
`MultiScalarField3D<int>` node mask holds three values: fluid (0), geometry
wall (1, from `isBoundary()` at set-up), clot wall (2, from
`PorousMedium::clotWall()`). `PorousMedium::updateWalls()` updates it for the
nodes it converts. The kernel sets the face velocity to zero on faces that
touch a node with mask 1 or 2, and the face diffusivity to zero for mask 1 and
to the factor of D9 times `D` for mask 2.

    class ScalarTransport {
    public:
      int addSpecies(const std::string & name, T diffusivity /* m^2/s */, T initialValue);
      void setBoundary(int species, int face /* 0..5 */, ScalarBoundary type, T value = 0);
      void addWallFlux(int species, plb::Box3D region, T flux);       // unit m^-2 s^-1
      plb::MultiScalarField3D<T> & concentration(int species);
      void step();                    // collective; one update of all species with dts = N_s dt
      void exchangeEnvelopes(int species);
      T totalAmount(int species);     // collective, sum over fluid and clot-wall nodes times dx^3
      void checkpoint();  static void restore(HemoCellFields &);
    };

The update loops are plain per-block loops over the bulk (a `HemoCellFunctional`
with `modif::nothing`; envelope exchange is explicit with
`duplicateOverlaps(field, modif::staticVariables)`, which is order-independent
for these fields because every envelope node has one owner).

### 4.4 `hemo::ReactionSystem` and `hemo::ReactionRule` (`helper/reactionRules.*`)

    struct ReactionState {            // one node, SI units
      const T * c;                    // [nSpecies]
      const T * phi;                  // [nComponents]
      T shearRate;                    // [1/s]
      T dt;                           // reaction step [s]
    };
    struct ReactionRates {            // rules add to these
      T * dcdt;                       // [nSpecies], unit/s
      T * dphidt;                     // [nComponents], 1/s
    };
    class ReactionRule {
    public:
      virtual ~ReactionRule() {}
      virtual void evaluate(const ReactionState & in, ReactionRates & out) const = 0;
      virtual std::string name() const = 0;
      virtual bool needsShearRate() const { return false; }
    };

`evaluate()` is called per node: no allocation, no logging, no look-up by name
(rules resolve species and component indices once, in their constructor).
`ReactionSystem::addRule()` takes ownership. The shear rate is computed only if
a rule asks for it. Built-in rules are created from XML by a factory; the case
calls `hemocell.reactions().addRule(new MyRule(...))` for its own.

### 4.5 Configuration

Proposed XML (new sections; absent sections mean the feature is off and nothing
is allocated):

    <porousMedium>
      <wallResistance> 1e3 </wallResistance>       <!-- alpha_wall, dimensionless (D8) -->
      <wallUpdateEvery> 20 </wallUpdateEvery>      <!-- iterations, multiple of the particle velocity update interval -->
      <protectCells> 1 </protectCells>             <!-- do not convert a node that holds a vertex -->
      <component>
        <name> platelet </name>
        <permeabilityLaw> kozenyCarman </permeabilityLaw>
        <diameter> 2.0e-6 </diameter>              <!-- m -->
        <bindsPlatelets> 1 </bindsPlatelets>
      </component>
      <component>
        <name> fibrin </name>
        <permeabilityLaw> davies </permeabilityLaw>
        <fibreRadius> 5.0e-8 </fibreRadius>        <!-- m -->
        <bindsPlatelets> 0 </bindsPlatelets>
      </component>
      <bindingSolidFraction> 0.3 </bindingSolidFraction>
      <repulsionSolidFraction> 0.5 </repulsionSolidFraction>
    </porousMedium>

    <scalars>
      <stepEvery> 4 </stepEvery>                    <!-- N_s, fluid steps -->
      <limiter> vanLeer </limiter>
      <cflMax> 0.3 </cflMax>
      <species>
        <name> thrombin </name>
        <diffusivity> 5.0e-11 </diffusivity>        <!-- m^2/s -->
        <clotWallDiffusivityFactor> 1 </clotWallDiffusivityFactor>   <!-- 0 makes clot walls impermeable (D9) -->
        <initial> 0 </initial>
        <boundary> <xNeg> dirichlet 0 </xNeg> <xPos> zeroGradient </xPos> </boundary>
      </species>
    </scalars>

    <reactions>
      <stepEvery> 20 </stepEvery>                   <!-- N_r, fluid steps -->
      <rule>
        <type> deposition </type>
        <component> fibrin </component>
        <species> thrombin </species>
        <rate> 50 </rate>                           <!-- 1/s -->
        <gate> step </gate> <threshold> 1e-6 </threshold>   <!-- mol/m^3 -->
        <shearMax> 1000 </shearMax>                 <!-- 1/s -->
        <phiMax> 0.2 </phiMax>
      </rule>
    </reactions>

The numbers above are placeholders for the format, not recommended values. The
reader in `config/config.*` iterates over repeated child elements; check how the
`Edge` list of `InnerEdges` is read and reuse that. Binding settings stay in the
cell XML (`<MaterialModel>`): `bindingDistance` [um], `shearThreshold` [1/s],
optional `shearMax` [1/s], `boundSolidFraction`. Keep reading the old
`distanceThreshold` (lattice units) with a deprecation warning.

All values are validated once at set-up (`fatal()` naming the element):
positive diffusivities and rates, known names, `phi` bounds, consistent
periodicity, intervals `>= 1`, `0 < wallResistance <=` the validated bound
(D8). At start-up the log states `alpha_wall` together with the permeability it
corresponds to (`nu dt / alpha_wall` in m^2) and, per component, the solid
fraction at which that component alone reaches it, so the user sees what will
turn into walls.

### 4.6 Output, checkpoints, restarts

* Fluid output constants (new numbers, do not renumber): `OUTPUT_SOLID_FRACTION`
  (one dataset `SolidFraction_<component>` per component, like
  `OUTPUT_CELL_DENSITY`), `OUTPUT_PERMEABILITY` (`K` in m^2 in SI output,
  `alpha` in lattice-unit output), `OUTPUT_SCALARS` (one dataset
  `Scalar_<species>` per species, in the Fluid files so they share the grid),
  `OUTPUT_CLOT_WALL` (1 in clot walls; `OUTPUT_BOUNDARY` keeps showing every
  wall). `OUTPUT_PERMEABILITY` is computed from `phi`, not read back from the
  lattice; in clot walls it is 0 (SI) or `alpha_max` (lattice units).
  `OUTPUT_BINDING_SITES` keeps its meaning. `OUTPUT_VELOCITY` needs no change:
  the dynamics returns the true velocity.
* Update `io/xdmfOutput.*` only if the attribute records need it (datasets
  written through `WriteFluidField` are listed automatically),
  `scripts/HDF5toXDMF.py`, and `doc/user_guide/reference/outputs.rst`.
* Checkpoint files in `<checkpointDirectory>`: `porous_<component>.{dat,plb}`,
  `bindingSites.{dat,plb}` (existing name), `scalar_<species>.{dat,plb}`;
  previous copies renamed to `.old` on rank 0 before saving, all before
  `checkpoint.xml`. `checkpoint.xml` `<General>` gains the component and species
  lists (checked at restart, fatal on mismatch), the bound-cell counters, and a
  format number, and `wallResistance` (a different value at restart is a
  warning: existing clot walls stay). `alpha` is recomputed from `phi` after
  loading, so it does not depend on what the lattice file holds. Clot walls
  come back twice: as bounce-back dynamics in the lattice file (as the legacy
  code relies on) and as `clotWall.{dat,plb}`; `restore()` checks that they
  agree (4.2).
* Restart on a different rank count must work; restart must be bit-for-bit on
  the same rank count.
* Load balancing: extend `loadBalancingSupported()` to refuse when the porous
  medium or scalars are enabled, and say so in `known_limitations.rst`.
* Pre-inlet: all new fields live on the domain lattice only. Every collective
  call is made on all ranks, pre-inlet ranks included (they own no blocks of
  these fields).
* `resetGlobalState()`: reset any new entry of `hemo::global`.

### 4.7 Transient binding (WP6)

* **Anchor list.** `std::vector<plb::Dot3D> anchorNodes` per
  `HemoCellParticleField`, next to `boundaryParticles`, block-local
  coordinates, bulk and envelope, sorted by node. Built by a function that
  `PorousMedium::refreshResistance()` and `updateWalls()` call with the other
  caches. A parallel vector holds the active flag, refreshed at each adhesion
  update.
* **Force.** `HemoCellParticleField::applyAdhesionForce()` and
  `HemoCellFields::applyAdhesionForce()` (a `HemoCellFunctional`), modelled on
  `applyBoundaryRepulsionForce()` (`core/hemoCellParticleField.cpp:809-838`):
  loop over the anchor list, visit the particle grid of the 3x3x3 nodes around
  each anchor, skip particles of cell types without tethering, add the two
  terms of 3.6 to `sv.force_repulsion`. The pair law is a free function
  declared in a header (compare `repulseParticlePair()` in
  `core/particleGrid.cpp`), so it can be unit-tested and later replaced (D10). Called in `iterate()` after the repulsion calls and before
  `capRepulsionForce()`, on repulsion update steps.
* **Arrest.** A function of `HemoCellParticleField` called from the binding
  step for the mode `tethered`, working on `get_particles_per_cell()` (complete
  cells only), writing to owned particles only. It needs, per vertex, whether
  the vertex is within range of an active anchor; compute that in the same call
  from the anchor list and the particle grid, not from a value stored by the
  force step.
* **Particle record.** `uint16_t arrestCount` after `solidify`
  (`core/hemoCellParticle.h`), zeroed in the constructors. A unit test asserts
  that `sizeof(serializeValues_t)` is unchanged; if it is not, stop and report
  (D11). Padding bytes are zero in existing checkpoints, so they read as a
  counter of 0.
* **Configuration**, per cell type inside `<MaterialModel>` of the cell XML,
  next to the settings of the `direct` mode (4.5), read once into
  `HemoCellField`, converted there to lattice units:

      <Binding>
        <mode> tethered </mode>                      <!-- direct (default) or tethered -->
        <tetherForce> 20 </tetherForce>              <!-- k_t, pN -->
        <tetherRange> 0.4 </tetherRange>             <!-- r_t, um; at most dx -->
        <tetherFriction> 0 </tetherFriction>         <!-- gamma_t, pN s/um; 0 switches the term off -->
        <tetherShearMin> 0 </tetherShearMin>         <!-- 1/s -->
        <tetherShearMax> 1e9 </tetherShearMax>       <!-- 1/s -->
        <minTetheredVertices> 3 </minTetheredVertices>
        <arrestSpeed> 1e-4 </arrestSpeed>            <!-- m/s -->
        <arrestTime> 1e-4 </arrestTime>              <!-- s -->
      </Binding>

  The numbers are placeholders for the format. Set-up checks (fatal): wall
  repulsion enabled, `0 < tetherRange <= dx`, `tetherForce >= 0`, the friction
  bound of 3.6, `arrestTime / (N_b dt)` below the range of the counter. The log
  states `k_t`, `r_t` and `gamma_t` in lattice units, the depth of the
  attraction well (`2/3 k_t r_t`) in units of `kBT`, and the friction bound.
* **Output.** No new field by default: adhesion shows in
  `OUTPUT_FORCE_REPULSION` (document this). Add `OUTPUT_ARREST_COUNT` as a
  per-vertex cell output, and the number of tethered cells per type as a column
  of `bound_cells.csv`.
* **Parallel.** The force is recomputed by every block from its own particles
  and its own anchor list, like repulsion, so it needs no communication beyond
  the existing envelope exchanges. The arrest evaluation follows the pattern of
  3.4: identical decision on every block, written by owners, spread by the next
  `syncEnvelopes()`.

---

## 5. Work packages

Order and dependencies:

    WP0 --> WP1 --> WP2 --+
       \--> WP3 ----------+--> WP4 --> WP5 --> WP6

WP1 and WP3 touch different files and can run in parallel after WP0. WP2 needs
WP1. WP4 needs WP1 and WP3. WP6 is scheduled last; technically it needs only
WP2. Each package ends with: unit tests green on 1 and 2
ranks (`ctest -L unit`), validation where physics changed, the `sanitize`
preset for the unit tests, `scripts/ci/smoke_drivers.py`, changed lines
formatted with `git clang-format`, user guide, `CHANGELOG.md` and `AGENTS.md`
updated. Report exactly what was built and run.

### WP0 Groundwork and decision gate

Goal: prove D1 on the real code base and put shared pieces in place.

Tasks

1. Add the descriptor (4.1) and select it as `DESCRIPTOR`. Build `hemocell`,
   `hemocell_test`, `hemocell_examples`, `hemocell_cases`. Fix what assumes
   three external scalars.
2. Add `GuoDarcyBGKdynamics` with `alpha` read from the external scalar, and
   the factory `hemo::newFluidDynamics(cfg)`.
3. Regression: run the pipeflow integration set-up (`tests/integration/pipeflow_setup.h`)
   with the plain Guo dynamics and with the new dynamics at `alpha = 0`;
   particle positions, velocities and the populations must be identical bit for
   bit after a few hundred steps on 1 and 2 ranks. Add this as a test.
4. Measure the cost of the wider cell and of the new dynamics on
   `tests/performance` (or pipeflow) with the same ranks and workload, using
   `hemo::global.statistics`. **Gate:** report both numbers to the maintainer.
   If the wider cell costs more than about 2 % of total run time, switch to
   D1(b): select the descriptor through a feature macro in
   `config/constant_defaults.h` and a variant `hemocell_porous`, add the macro
   to `HEMOCELL_BUILD_CONFIGURATION`, and build a second test runner for it.
5. Lattice layout in the checkpoint fingerprint, fatal on mismatch (4.1).
6. D3: verify with `static_assert`/a unit test that `sizeof(serializeValues_t)`
   is unchanged when `solidify` is always present; if so, remove the
   `SOLIDIFY_MECHANICS` conditionals, keep `hemocell_solidify_mechanics` as an
   alias of `hemocell` so existing case `CMakeLists.txt` keep working, and
   adjust the fingerprint string and `HEMOCELL_BUILD_CONFIGURATION`. If the
   size changes, stop and report.
7. `hemo::shearRate()` (3.3) in `helper/`, with unit tests in Couette and
   Poiseuille flow (analytic `gamma_dot`, 1 % in the bulk away from walls).
   Do not change the binding criterion yet.

Acceptance: all existing tests pass; `alpha = 0` regression is bit-identical;
performance numbers reported; old-format checkpoints are refused with a clear
message.

### WP1 Porous medium

Goal: a fluid that feels a prescribed solid-fraction field.

Tasks

1. `PorousMedium` (4.2) with components, permeability laws (unit-tested as pure
   functions against hand-computed values), `refreshResistance()`,
   `exchangeEnvelopes()`.
2. The wall switch (3.1.1): `clotWall` field, `updateWalls()` with the
   decision in the bulk, envelope exchange, bounce-back on bulk and envelope,
   and the cache updates of rule 7 (the scalar mask hook is a no-op until WP3).
   The `protectCells` guard is written here and tested with cells in WP2.
   Wire `updateWalls()` into `iterate()` and the interval into `sanityCheck()`.
3. Run V9 and record the validated bound of `alpha_wall` as a named constant
   next to the dynamics, with the test that established it. If the scheme does
   not pass up to the default `1e3`, stop and report the largest value that
   does (D8).
4. XML section `<porousMedium>`, validation, `HemoCell` accessor
   (`hemocell.porous()`, fatal if not enabled), creation in
   `initializeCellfield()`.
5. Output constants and writers (4.6); checkpoint and restore.
6. A setter for tests and cases: fill a component inside a box or from a
   functional of the node position.
7. Validation tests (section 6: V1 to V9).
8. Docs: a concept page (model, equations, units, limits, the table of 3.1.1
   for choosing `wallResistance`), `xml_files.rst`, `outputs.rst`,
   `hemocell_api.rst`, `known_limitations.rst`.

Acceptance: V1 to V9 pass on 1 and 2 ranks; restart is bit-identical, with
clot walls and with pending nodes; a run with the porous medium enabled and all
`phi = 0` is bit-identical to one without; the validated bound of `alpha_wall`
is reported.

### WP2 Binding deposits solid fraction

Goal: feature 1 of the request, correct on several ranks.

Tasks

1. Rewrite the binding step as in 3.4: owner-only flagging, deposit in the
   bulk only, removal on all blocks, envelope exchange, refresh. Keep the
   `CellMechanics` hook but give it the porous medium instead of the CEPAC
   block; rename it (`bindingMechanics`) and keep `solidifyMechanics` as a
   deprecated forwarder if that is cheap.
2. Read thresholds once at set-up into `HemoCellField`; units as in 4.5; use
   `hemo::shearRate()`. This changes the meaning of `shearThreshold` (defect 3):
   state it in `CHANGELOG.md` under "fixes that change simulation results" and
   update the shipped `PLT.xml` files.
3. Legacy behaviour through the `impermeable` law and the wall switch (D4); no
   second code path. Remove the direct `defineDynamics(..., BounceBack)` from
   the mechanics model.
4. Dense drag nodes and clot walls in the wall-repulsion list (3.4); the
   `protectCells` guard exercised with cells; `sanityCheck()` for the intervals
   (defect 7) and a warning when the porous medium runs without wall repulsion.
5. Counters, CSV line per output step (`bound_cells.csv`: iteration, time,
   bound cells per type, solid volume per component), checkpointed.
6. Tests B1 to B6 (section 6). Migrate `cases/solidify_example`,
   `cases/fc_solidified_stenosis`, `cases/CEPAC` (binding part) and the smoke
   script entries.
7. Close the `TODO.md` item "Solidification communication" or restate what
   remains.

Acceptance: B1 to B6 pass; 1-rank and 2-rank runs of B2 give the same solid
fraction field, clot walls and bound-cell count; the two cases run under the
smoke script.

### WP3 Scalar transport

Goal: feature 2 of the request. Independent of WP1/WP2 apart from WP0.

Tasks

1. `ScalarTransport` (4.3): fields, velocity field, the three-valued node mask
   and its update hook for `PorousMedium::updateWalls()` (the hook is tested
   here by setting mask values directly, so WP3 does not depend on WP1), the
   update kernel with limiter options, boundary conditions, sources.
2. XML section `<scalars>`, validation, stability checks, `hemocell.scalars()`.
3. Output (`OUTPUT_SCALARS`), checkpoint, restore, total-amount logging.
4. Tests S1 to S9 (section 6).
5. Port `cases/CEPAC` to the new solver; then remove the CEPAC lattice:
   `enableCEPACfield`, `createCEPACfield`, `CEPACfield`, `CEPAClattice`,
   `CEPAC_DESCRIPTOR`, `param::tau_CEPAC`, `setCEPACOutputs`,
   `writeCEPACField_HDF5`, the checkpoint branch, the load-balancer condition,
   and their documentation. `<enableCEPACfield>1` in a configuration becomes a
   fatal error that names the replacement. List the removal under "Structure"
   in `CHANGELOG.md`.
6. Docs: concept page (equation, scheme, stability limits, what is not
   resolved), `xml_files.rst`, `outputs.rst`, API reference.

Acceptance: S1 to S9 pass on 1 and 2 ranks; cost per species per update
reported from the profiler; no CEPAC symbol left outside `legacy/` and the
changelog.

### WP4 Reaction rules

Goal: feature 3 of the request.

Tasks

1. `ReactionRule`, `ReactionSystem`, built-in rules and their XML factory (4.4, 4.5).
2. The reaction step in `iterate()` (3.5): gather inputs per node (fluid and
   clot-wall nodes, shear rate 0 in the latter), call rules, update, clamp,
   exchange, refresh `alpha` and the binding and repulsion lists. The step does
   not convert nodes; nodes it pushes past `alpha_wall` are pending until the
   next wall update.
3. Rule sources feed `ScalarTransport` either directly in the reaction step
   (operator splitting, default) or as a source field; choose one and document it.
4. Tests R1 to R6 (section 6).
5. Docs: a developer recipe "Adding a reaction rule" in `developer/extending.rst`;
   user page with every built-in rule, its equation, parameters and units.

Acceptance: R1 to R6 pass; a custom rule written in a test compiles against the
public headers only.

### WP5 Integrated example and clean-up

1. `cases/thrombus_example`: a small channel in shear flow with an injury patch
   on the wall: thrombin wall flux on the patch, platelets binding to it,
   bound platelets releasing an agonist, fibrin deposition where thrombin
   exceeds a threshold at low shear, permeability falling accordingly. Small
   enough for the smoke script. Register it in `cases/CMakeLists.txt`, the root
   `.gitignore` and `scripts/ci/smoke_drivers.py`. State in its README and in
   the user guide that the parameters are illustrative, not calibrated.
2. A user-guide page for the case with what to look at in ParaView.
3. `known_limitations.rst`: scalar passes through membranes; concentration
   boundary layers at reactive walls are thinner than `dx`, so wall mass
   transfer is grid dependent; clot surface is a staircase; material less
   permeable than `nu dt / alpha_wall` carries no flow, and a clot wall does
   not turn back into fluid; the flow rate jumps slightly when a layer converts
   (measured in V4); BGK at the working `tau` (error measured in V2, V4); no
   load balancing; not combined with interior viscosity; bound platelets are
   removed, so clot mechanics (deformation, embolisation) are not represented.
4. `AGENTS.md`: source table, variants, `iterate()` order, checkpoint contents,
   verification matrix.
5. Run the full matrix: unit and validation on 1 and 2 ranks, `sanitize` and
   `debug` presets for the unit tests, smoke run with `--slow`.

### WP6 Transient binding

Goal: tethering, slowed motion along the surface and arrest before the deposit
(3.6, 4.7). Needs WP2.

Tasks

1. Anchor node list and its maintenance with the other caches; active flags
   from `hemo::shearRate()`.
2. The pair law as a free function with unit tests; `applyAdhesionForce()`;
   the call in `iterate()`; `sanityCheck()` entries (wall repulsion enabled,
   intervals).
3. `arrestCount` in the particle record with the size test (D11); the arrest
   evaluation in the binding step; binding mode per cell type, `direct`
   unchanged.
4. Cell XML `<Binding>` block, validation, unit conversion, start-up log (4.7).
   Add the block to the platelet template in `tools/cell_templates`.
5. `OUTPUT_ARREST_COUNT`, the tethered-cell column of `bound_cells.csv`.
6. Tests T1 to T8 (section 6).
7. Extend `cases/thrombus_example` with a tethered platelet type, or add
   `cases/platelet_tethering` (one platelet in shear flow over a binding
   patch) if that keeps the example small; register it as in WP5.
8. Docs: model page (anchor nodes, both force terms, arrest criterion, what the
   model does not represent), `xml_files.rst`, `outputs.rst`,
   `known_limitations.rst`, `CHANGELOG.md`, `AGENTS.md`.

Acceptance: T1 to T8 pass on 1 and 2 ranks; with no cell type in the mode
`tethered` a run is bit-identical to one before WP6; restart during tethering
is bit-identical.

---

## 6. Tests and validation

Fluid-only tests need no cells; build a small box or channel directly (see
`tests/integration/` for set-up patterns, `hemo_test::testConfig()` for the
configuration, `HemoCell::MPIHandle::External`). All formulas are in lattice
units unless a unit is given. Tolerances below are starting points derived from
the discretisation; tighten them where the result allows, and do not widen them
to get a pass.

Porous medium

| ID | Case | Check |
|---|---|---|
| V1 | Periodic box, uniform `alpha`, body acceleration `g` | Steady `u = g / alpha`, relative error < 1e-10, for `alpha` = 1e-3, 1, 1e3 and two values of `tau`; before steady state, the per-step factor of 3.1 |
| V2 | Channel with bounce-back walls at `z = +-h`, uniform `alpha`, body acceleration `g` | `u(z) = (g/alpha) [1 - cosh(z/delta) / cosh(h/delta)]`, `delta = sqrt(nu_lbm/alpha)`; with `delta >= 3` nodes, L2 error < 2 %; record the error at the production `tau` (about 1.8) |
| V3 | Channel half free, half porous | Continuity of velocity; flow rate against the analytic two-layer solution, 3 % |
| V4 | Poiseuille flow between two slabs, three ways: drag nodes just below `alpha_wall` (conversion disabled), the same slabs converted to clot walls, geometry bounce-back walls | Clot walls and geometry walls give identical flow, bit for bit; velocity inside the drag slabs < `2 / alpha_wall` of the peak; report the effective wall position of the drag slabs and the jump in flow rate at conversion, which must lie within what a half-node wall shift explains |
| V5 | `alpha = 0` everywhere, pipeflow with cells | Bit-identical to the plain Guo dynamics (WP0) |
| V6 | Step change of `alpha` in a running flow, from 1e-3 up to `alpha_wall`, and pending nodes (conversion held back) up to `alpha_max` | No NaN; `|u|` decreases monotonically in the porous region; density stays within a stated bound of its mean |
| V7 | Porous block across a block boundary, 1 and 2 ranks; checkpoint and restart | Identical fields; restart bit-identical |
| V8 | Channel flow in which `phi` of a block of nodes is raised in steps past the threshold; the block straddles a block boundary and a periodic boundary; 1 and 2 ranks; checkpoint with pending nodes and with clot walls | Nodes convert at the first wall update after `alpha >= alpha_wall`, not before; the lattice has bounce-back dynamics on exactly the flagged nodes, in bulk and envelope of every block; populations identical on 1 and 2 ranks; density disturbance at conversion below a stated bound; restart bit-identical; every cache of rule 7 equals a rebuild from scratch; a later decrease of `phi` leaves the wall in place |
| V9 | V1, V2 and V6 repeated for `alpha` = 2, 10, 1e2, 1e3, 1e4 at `tau` = 1 and the production `tau` | The largest `alpha` for which all three pass is the validated bound of `alpha_wall` (D8); record the Darcy velocity error and the decay of the oscillating mode at each value |

Scalars

| ID | Case | Check |
|---|---|---|
| S1 | Uniform `c` in pipeflow with cells | Stays uniform to 1e-12 |
| S2 | Gaussian, no flow, periodic box | Variance grows as `sigma0^2 + 2 D t`, 2 % |
| S3 | Gaussian (`sigma` = 4 nodes) and step, uniform flow, one pass through a periodic box | Peak retained > 90 % (Gaussian); `min >= -1e-12`, `max <= c0 + 1e-12` (step, TVD limiters); L1 error at least 5 times smaller than first-order upwind |
| S4 | Wall kept at fixed `c` in shear flow (Leveque problem) | Local wall flux against the Leveque solution over the resolved range; report where the boundary layer is thinner than two nodes |
| S5 | Closed box with walls, with and without flow, and a periodic box with cells | Total amount constant to round-off without flow; drift with flow reported and below a stated bound |
| S6 | Any of the above on 1 and 2 ranks | Identical fields (fluid-only cases bit for bit) |
| S7 | Checkpoint and restart, same and different rank count | Bit-identical on the same rank count |
| S8 | `N_s` = 1, 2, 4 | Results converge; stability check fires when `CFL_max` is exceeded |
| S9 | Two reservoirs at fixed `c` separated by a slab, no flow: slab of clot-wall nodes (factor 1 and 0.5), slab of geometry-wall nodes; then a fluid node that converts while it holds scalar | Steady flux `D_slab * delta_c / L` to 1 % through the clot wall; exactly zero through the geometry wall and with factor 0; the converting node keeps its amount and the total is unchanged by the conversion |

Reactions

| ID | Case | Check |
|---|---|---|
| R1 | Uniform fields, no flow, `deposition` without gate | `phi_n = phi_max (1 - (1 - k dtr)^n)` to 1e-12; analytic `1 - exp(-k t)` within the Euler error |
| R2 | Concentration below and above the step threshold | Exactly no change below |
| R3 | Couette flow with known `gamma_dot`; shear window | Rule active only inside the window; `hemo::shearRate()` within 1 % of the analytic value in the bulk |
| R4 | Channel flow, thrombin source upstream, `deposition` of fibrin; once with `phi_max` below and once above the solid fraction that reaches `alpha_wall` | Fibrin appears only where `c >= c_thr`; below the threshold the flow rate falls towards the Darcy value of V1/V2 for the final `phi`; above it the nodes become clot walls at a wall update, the rule keeps acting in them with shear rate 0, and the result is the same on 1 and 2 ranks |
| R5 | Sink larger than the available amount; `phi` reaching `phi_max` | No negative values, clamps counted and logged |
| R6 | 1 and 2 ranks; checkpoint and restart | Identical; bit-identical restart |

Binding

| ID | Case | Check |
|---|---|---|
| B1 | One platelet in shear flow next to a wall binding site, shear inside the window; once with a permeability law that stays below `alpha_wall`, once with `impermeable` | It binds; its particles are gone on all blocks; deposited nodes equal the inner nodes; `phi` and the binding flags are set there. Below the threshold the nodes stay drag nodes with the expected `alpha`; with `impermeable` they are clot walls after the same step |
| B2 | Same with the platelet straddling a block boundary, 1 and 2 ranks (and a periodic boundary) | Same solid-fraction field and counters as on 1 rank; no node written twice |
| B3 | Second platelet arriving at the first deposit | Binds to deposited nodes |
| B4 | Checkpoint during aggregation, restart | Bit-identical continuation |
| B5 | RBC passing a dense deposit with repulsion enabled | The RBC is not deleted; no vertex ends inside nodes with `phi >= phi_repulsion` beyond the repulsion cutoff |
| B6 | An RBC resting on nodes whose `phi` is raised past the threshold, with and without `protectCells`; 1 and 2 ranks with the RBC across the block boundary | With the guard: occupied nodes stay pending and are counted, the RBC survives, and the nodes convert at the first wall update after it has left. Without: the RBC is deleted and the deletion is logged. Same result on 1 and 2 ranks |

Transient binding

| ID | Case | Check |
|---|---|---|
| T1 | Pair law as a function | Zero at `d = 0` and `d = r_t`, `k_t` at `r_t / 2`, direction towards the anchor; friction weight 1 at the anchor and 0 at `r_t`; nothing beyond `r_t` |
| T2 | One platelet at rest above a binding wall, no flow | It settles with its near side where attraction and wall repulsion balance (computed for the test's parameters); no vertex is deleted; the anchor list equals the non-wall face neighbours of the binding nodes, on every block including envelopes |
| T3 | The same platelet in shear flow, `gamma_t = 0`, shear rate varied | Mean speed along the wall is below that of an identical run without adhesion; below some shear rate it is pinned; report that critical shear rate and its dependence on the flow direction relative to the lattice axes |
| T4 | As T3 above the critical shear, `gamma_t` varied | Speed decreases monotonically with `gamma_t`; the set of tethered vertices changes while the platelet moves and its orientation changes (rolling), which distinguishes it from rigid sliding; the stability check fires above the bound |
| T5 | Arrest | With speed below `arrestSpeed` for `arrestTime` the platelet is deposited exactly as in B1, at the step the counter predicts; a platelet that stays faster, or that loses its tether, is not, and its counter returns to 0 |
| T6 | Tethered platelet crossing a block boundary and a periodic boundary, 1 and 2 ranks | No jump in force or position at the crossing; same arrest step; trajectories agree within the tolerance that cell runs on different decompositions have without adhesion |
| T7 | RBCs only, adhesion configured for a platelet type that is absent; then RBCs passing a tethered platelet | First run bit-identical to one without adhesion; in the second no RBC vertex receives adhesion force |
| T8 | Checkpoint while tethered with a non-zero counter, restart | Bit-identical continuation; a checkpoint written before WP6 loads with counters 0 |

Regression and cost

* A run with every new feature compiled in but disabled is bit-identical to the
  run before this work (pipeflow, 1 and 2 ranks, checkpoint files compared with
  `cmp` where the format is unchanged, otherwise particle data).
* `Output.DoesNotChangeTrajectory` extended to runs with the porous medium,
  scalars and reactions enabled.
* Profiler numbers for: wider cell, new dynamics, scalar step per species,
  reaction step per rule, binding step. Same case, ranks and workload.

---

## 7. Rules for the implementing agents

* One work package per branch or pull request; do not start a package whose
  prerequisites are not merged.
* Lattice dynamics are changed in one place only, `PorousMedium::updateWalls()`,
  on bulk and envelope, at steps with synchronised particle envelopes.
* Stop and report instead of working around: the WP0 gate, V9 not confirming
  the default `alpha_wall`, a change of
  `sizeof(serializeValues_t)`, a validation that fails by more than its
  tolerance, any need to edit `palabos/` (integration changes go to
  `patch/palabos.patch`; none is expected).
* Library code is C++11. No new dependency.
* Block-level loops: no allocation, no map look-ups, no XML reads, no logging.
  Hoist `param::` values. Skip nodes with `alpha = 0`/`phi = 0` where that
  avoids work, without branches in the collision itself.
* Every new collective call is made by all ranks in the same order, pre-inlet
  ranks included. Never call `MultiScalarField3D::get()` in per-rank code.
* Envelope exchange of the new fields is explicit, right after the field
  changes. Nothing may depend on message arrival order.
* Units at interfaces: micrometres for distances in cell XML, metres and SI
  elsewhere, iterations for intervals. Convert once, where the value enters.
* Output and statistics do not change state.
* Public API and XML tags get Doxygen comments and user-guide entries that
  state units, valid ranges and whether the call is collective. Do not document
  validation that the code does not perform.
* Shipped parameter values that are not taken from a cited source are labelled
  as illustrative.

---

## 8. Risks and open questions

| Risk | Consequence | Mitigation |
|---|---|---|
| Custom descriptor interacts badly with a Palabos template (specialisations, boundary conditions, pre-inlet, voxeliser) | WP0 stalls | WP0 builds every driver first; fallback D1(b) or per-node dynamics |
| BGK accuracy of the Brinkman term at `tau` about 1.8 | Effective permeability and wall position depend on `tau` | V2 and V4 measure it; a TRT collision is a later, separate change |
| Flux-limited unsplit scheme is not strictly monotone in 3D | Small over- or undershoots | S3 bounds them; `CFL_max` guard; dimensional splitting as an option if the bounds fail |
| Concentration boundary layers thinner than `dx` at reactive walls | Wall fluxes and thresholds are grid dependent | Document; S4 shows the resolved range; subgrid wall model is out of scope |
| Scalars cross cell membranes | Agonists diffuse through RBCs | Document; a membrane-aware scheme is out of scope |
| Thin platelets have few inner nodes | Deposited volume differs from platelet volume | Fallback node (3.4); log deposited volume against bound-cell volume |
| Time-scale gap between flow and biochemistry | Cases need scaled rates | Rates are explicit inputs; example states its scaling |
| Removal of the CEPAC lattice and of the solidify variant breaks private cases | Users must port | Fatal errors that name the replacement; alias target; changelog |
| `alpha_wall` set low (near 2) | Most clot material becomes impermeable; feature 1 behaves like the present bounce-back model | Default `1e3`; start-up log states the permeability and solid fractions that the threshold corresponds to (4.5) |
| Conversion changes the effective wall position by about half a node | Small jump in flow rate and a density disturbance when a layer converts | V4 and V8 bound both; conversions are spread over wall updates as the clot grows |
| A node converts while a cell is in or next to it | Cell loses vertices and is deleted | Repulsion from dense nodes before conversion; `protectCells` guard; B6 |
| Pending nodes never convert because a trapped cell stays in them | A porous pocket with a cell inside a wall region | Pending count is logged; cases can switch the guard off |
| Clot walls are one-way | Lysis or detachment cannot reopen a node | Stated as a limitation; see open question 1 |
| Interior viscosity resets node dynamics | A clot wall could be turned back into fluid | Combination refused at set-up until tested (4.2) |
| Anchors form a grid with spacing `dx` | Pinning and sliding depend on the direction of the flow relative to the lattice and on the staircase of curved surfaces | T3 reports the dependence; friction, not the corrugation, is the intended source of slowing, so keep `k_t` moderate |
| Attraction pulls vertices close to the wall | A vertex crosses into a wall node and its cell is deleted | Anchors lie outside the surface; wall repulsion is mandatory with adhesion; T2 |
| `FORCE_LIMIT` caps repulsion and adhesion together | Strong adhesion is silently limited | Start-up log compares `k_t` with the limit; warning when `k_t` exceeds half of it |
| Friction with a lagged velocity | Oscillation or instability for large `gamma_t` or long update intervals | Bound checked at set-up (3.6); T4 |
| Near-wall lubrication is not resolved | Speed along the wall depends on the stand-off, not on a physical film | Documented; calibrate `gamma_t` against measured translocation speeds |
| Translocation takes seconds, runs cover milliseconds | Arrest thresholds have to be scaled | Thresholds are explicit inputs; example states its scaling |

Open questions for the maintainer, none blocking WP0:

1. Should bound platelets be able to detach (solid fraction returning to
   cells, or simply being removed under high shear)? The design allows `phi` to
   decrease by a rule but does not re-create cells, and a clot wall does not
   revert. Reverting needs a rule for the populations of the reopened node
   (for example equilibrium at rest with the mean density of its fluid
   neighbours) and a `clotWall` flag that can be cleared; neither is planned.
2. Should the scalar diffusivity depend on the solid fraction from the start?
   The kernel takes a per-face diffusivity, so it can be added without a
   redesign.
3. Is a Forchheimer (quadratic) drag term wanted? The implicit velocity formula
   has a closed form for it as well.
4. Per-cell state beyond the binding flag and the arrest counter (activation
   level, bound time) would change the particle record; it is not part of this
   plan.
5. Should the anchors sit exactly on the first non-wall node (3.6), or should
   the distance from the surface be a parameter? A parameter means anchors are
   points instead of nodes, which costs a separate list of positions but
   removes the tie between the stand-off and `dx`.
6. Should tethered platelets also attract each other before either is bound
   (aggregation in the flow)? The design tethers only to binding nodes.

---

## 9. Out of scope

Clot mechanics and embolisation, platelet activation state, discrete stochastic
bonds with force-dependent rupture (3.6), moving binding surfaces,
scalar-dependent cell mechanics, membrane-impermeable scalar transport, subgrid wall
mass-transfer models, TRT or MRT fluid collision, load balancing with the new
fields, calibrated biochemical parameter sets.
