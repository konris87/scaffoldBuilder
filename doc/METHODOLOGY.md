# ScaffoldBuilder — Algorithm & Methodology

**Living document.** Update it whenever the algorithm, a parameter, or a metric
changes. Each section states *what* is computed, *on which substrate* (mesh vs.
generation grid vs. resampled image), and any known bias.

Last updated: 2026-07-26

**Contents:** §1–7 algorithm & metrics · §8 verification/calibration tools ·
§9 empirical findings · §10 limitations · §11 changelog · §12 **experiment
design** (paper results structure).

---

## 1. Pipeline overview

```
container ──> seeds ──> scalar field ──> post-processing ──> iso-surface ──> metrics
 (SDF)      (Poisson)   (Voronoi +        (Taubin, clamp,    (Lewiner MC33)
                         anisotropy)       islands, seal)
```

Two distinct resolutions exist and must not be confused:

| name | symbol | meaning |
|---|---|---|
| **generation grid** | `stepX/Y/Z`, `blockDims` | where the scalar field is sampled and the mesh is extracted |
| **measurement voxel** | `voxelSize` | the target (µCT) resolution the structure is *resampled* to for voxel metrics |

**Rule:** every voxel-based metric must be measured at the *same* `voxelSize`,
matched to the reference µCT resolution.

---

## 2. Container

Three types (`box`, `cylinder`, `abstract mesh`), each exposing a signed
distance function (SDF). The mesh container's SDF uses the angle-weighted
pseudonormal with an AABB hierarchy for the nearest-triangle query.

The SDF is used to (a) confine seeds, (b) clamp the lattice to the container,
and (c) optionally drive graded thickness.

---

## 3. Seed generation

- **Random**: N uniform points in the container's bounding box.
- **Uniform Poisson**: minimum-distance constraint `r` (Bridson). Grid cell
  `cs = r/√3`; candidates drawn in the annulus `[r, 2r]`.
- **Varied Poisson**: the minimum distance is a function
  `R = ρ(Ψ(x); R_min, R_max)` of a user-chosen SDF `Ψ` (point / plane /
  container) and a univariate function `ρ` (linear, quadratic, smoothstep,
  constant).
- **Stochastic Poisson** (`stochasticRadius`): no distance field — each seed
  draws its own minimum distance from a **truncated normal** `N(μ, σ)` clamped
  to `[R_min, R_max]` (rejection-sampled to keep the Gaussian shape, not clamp
  mass to the bounds; `μ ≤ 0` defaults to the range midpoint). Variable-radius
  Bridson: the grid stays at `cs = R_min/√3` and each accepted seed stamps a
  neighbourhood of `⌈r/R_min⌉` cells so larger-`r` seeds are still tested. All
  draws come from the run's single seeded RNG stream, so a given `rngSeed`
  reproduces the whole cloud (positions **and** radii). Purpose: inject
  **intra-sample** spacing heterogeneity (a spread of cell sizes) for realism —
  e.g. toward a literature Tb.Sp *distribution*.

  ⚠️ **Input `σ` is not output Tb.Sp SD.** Two effects break the identity:
  radius→Tb.Sp is a monotone *transfer* (not 1:1), and the packing **regularises**
  the injected variance. Measured (4 mm box, `[0.55, 0.85]`): raising `σ`
  0.02→0.12 *lowered* mean Tb.Sp 0.85→0.78 (more small-`r` candidates pass the
  acceptance test → denser) while the intra-sample Tb.Sp SD barely moved
  (0.157→0.142). So `σ` is a **heterogeneity knob**, the **mean** still needs the
  usual radius/`μ` calibration, and the input→output transfer must be read off the
  profiler `--sweep-radius-std` map (§8.5), not assumed. Also note the spread this
  controls is **intra-specimen**, which is *not* the same SD as a literature
  inter-specimen `Tb.Sp ± SD` (§9, "three SDs").

Seed spacing is the dominant control on **Tb.Sp** (and hence cell size).

---

## 4. Anisotropy

### 4.1 Representation — covariance `C`, metric `M`

Each source stores a **covariance** (not a metric):

```
C = Rᵀ S² R          (R orthogonal, S = diag(sx, sy, sz))
```

`R` is built from a material direction (mapped to local X) plus a roll angle.
The **metric** used to measure distance is its inverse:

```
M = C⁻¹   ( = Rᵀ S⁻² R )
```

Stretch values are clamped (`max(1e-6, s²)`) so `C` can never be singular.

> **Naming:** `C` is what we *blend*; `M = C⁻¹` is what we *measure with*.
> The uniform path warps seeds by `S⁻¹R`, a square root of `M`
> (`‖S⁻¹Rδ‖² = δᵀMδ`) — the two paths are mathematically identical.

### 4.2 Blending — partition of unity

Each source contributes a **compactly supported** Gaussian weight:

```
w_m(x) = max(0, exp(−‖x − o_m‖² / 2σ_m²) − e^{−9/2})
```

The `e^{−9/2}` shift makes the weight vanish *continuously* at the `3σ_m`
support boundary; an unshifted hard cut-off would jump by ≈1.1% there and
imprint a spherical seam. `3σ` is the conventional truncation radius: ~1%
of peak discarded, ~1% shift distortion, and support cost grows as the cube.

The background acts as an ever-present source of fixed weight `w_b`:

```
C(x) = ( Σ_m w_m C_m + w_b C_bg ) / ( Σ_m w_m + w_b )
M(x) = C(x)⁻¹
```

Far from all sources this reduces *exactly* to `C_bg`; at a source centre the
background retains a residual `w_b/(1+w_b)`. `w_b` is clamped away from zero
(at `w_b = 0` the blend is undefined outside all supports).

### 4.3 Nearest-seed query

- **Uniform (`A = ∅`)**: seeds warped once by `S⁻¹R`, ordinary Euclidean kd-tree.
- **Varied (`A ≠ ∅`)**: unwarped kd-tree returns `K` Euclidean candidates
  (`K = min(max(⌈3κ²⌉, 12), |S|)`, `κ = s_max/s_min`), re-ranked by
  `d² = δᵀM(x)δ`. This is a heuristic, not a guarantee — the metric-nearest
  seed can in principle fall outside the Euclidean candidate set.

---

## 5. Scalar field

### 5.1 Domain and bands

The grid covers the container AABB. All *numerical* bands are expressed in
**grid spacings**, never absolute mm, so they track model scale and resolution:

| band | value | purpose |
|---|---|---|
| `surfaceMargin` | `3h` | shell of valid lattice kept outside the wall so clamp/smoothing/MC see a smooth gradient |
| `air_skip_level()` | `isoLevel + 2h` | value stamped into exterior voxels **and** the smoothing skip threshold (must be identical — the skip test is `>=`) |

where `h = max(stepX, stepY, stepZ)`.

### 5.2 The two geometric extremes

```
wallVal  = d₂ − d₁     zero on Voronoi FACES  -> thickened gives PLATES
strutVal = d₃ − d₁     zero on Voronoi EDGES  -> thickened gives RODS
```
Note `wallVal ≤ strutVal` always.

### 5.3 Openness and spread (per-face fenestration)

A Voronoi face is identified by its **nearest seed pair**. Hashing that pair
gives each face a fixed pseudo-random `faceP ∈ [0,1)` that every voxel near the
face agrees on — a coherent membrane, no storage, thread-safe, reproducible.

```
localTau = clamp( threshold + spread · (faceP − 0.5), 0, 1 )
value    = (1 − localTau)·wallVal  + localTau·strutVal
```

| parameter | meaning |
|---|---|
| **`threshold` (openness)** | **mean** openness. `0` → foam (plates on every face); `1` → lattice (bare rods). Drives **BV/TV**. |
| **`spread`** | per-face **diversity** of openness. `0` → every face identical (the plain uniform morph); `1` → faces range from near-full plate to bare rod (rod-and-plate coexistence). |

An intermediate `localTau` eats a plate back toward its bounding edges — i.e. a
**fenestrated** plate. Different faces are perforated to different degrees.

Properties:
- `E[localTau] = threshold` exactly (faceP is uniform), so `spread` reshapes the
  mixture without shifting mean openness.
- **The field is continuous across Voronoi edges** despite `localTau` jumping
  between faces: at an edge `d₁=d₂=d₃`, so `wallVal = strutVal = 0` and
  `value = 0` regardless of `localTau`. Only the gradient kinks; Taubin absorbs it.
- **Linear blend by default.** The base field is a linear blend, not a
  `min()`/union, so there is no crease to fillet and no `smin` is applied. A
  junction fillet is available as an **opt-in** (§5.3b); left off, rod–plate
  junctions stay sharp and the local-thickness metric does not amplify them via
  junction bleeding (see §7.2).

### 5.3b Optional junction smoothing — the `k` fillet

The linear blend of §5.3 is crease-free but its **rod–plate junctions are sharp**:
where a fenestrated plate meets the bare rod on a Voronoi edge the two surfaces
join at an angle. Real trabecular nodes are rounded, thickened blends. Turning on
`smoothJunctions` applies a **smooth-minimum fillet** (Inigo Quilez's polynomial
`smin`, `smin_gradient()`) between the linear-blend field and a *shifted strut
field*:

```
rodField = strutVal + (1 − localTau)·(strutVal − wallVal)
value    = smin( value, rodField, k )          // k = smoothK
```

`rodField ≥ value` everywhere, and their gap `2(1 − localTau)(strutVal − wallVal)`
→ 0 on Voronoi **edges** (there `wallVal = strutVal = 0`). So the fillet appears
**exactly at the junctions** and vanishes elsewhere — it does not move the plates
or the edges, only rounds the corner between them.

The polynomial `smin(a, b, k)` (value **and** analytic gradient, so the Newton
normalization of §5.4 stays valid):

```
if |a − b| ≥ k :   return the hard min (no blend)
else           :   h   = (k − |a − b|) / k          ∈ (0,1]
                   val = min(a,b) − h³·k / 6         (dips BELOW the min → fattens solid)
                   grad = blend(gradA, gradB; m = ½h²)
```

| parameter | meaning |
|---|---|
| **`k` (`smoothK`)** | **fillet radius**, in field units. Blending is active only where the two fields differ by less than `k`; `k → 0` reduces *exactly* to the linear blend (sharp), and larger `k` gives a **wider, rounder** node. Off by default (`smoothJunctions = false`). |

> ⚠️ The inline `@param k` comment on `smin_gradient()` ("higher = sharper") is
> **inverted**; the code behaves as tabulated above (higher `k` = smoother).

Because the fillet **fattens** the nodes, it lifts the *minimum* achievable
Tb.Th: the smin dips the field below `min(a,b)`, so more voxels fall under the
iso-level. That inflation is absorbed by `calibrate_thickness` (§8.2) — but if a
target Tb.Th sits below the fattened floor, the calibration cannot reach it, and
the joint solve **backs `k` off** until it can (§8.2). The extra Tb.Th is real
geometry, the same effect as junction bleeding (§7.2), which is why the fillet is
opt-in and kept small (`k ≈ 0.01` with typical thicknesses).

### 5.4 Gradient normalization

`value` scales with seed density and metric distortion, so `c_glob` is not a
length. Recover one with a single Newton step toward the surface:

```
D(x) ≈ v(x) / ‖∇v(x)‖          (Euclidean distance to the medial set)
Φ(x) = 2·v/‖∇v‖ − c_loc + c_glob
```

The solid is everything within `c_loc/2` of the medial set — a wall of full
width `c_loc` — so the factor of two converts the one-sided distance into a full
thickness, and the surface extracted at `c_glob` reproduces it.

The per-seed gradient of `d_k = (δᵀM(x)δ)^{1/2}` is

```
∇d_k = M(x)δ/d_k  +  (1/2d_k)·[ δᵀ ∂ᵢM(x) δ ]ᵢ
```

| regime | metric-derivative term | consequence |
|---|---|---|
| **uniform** (`A = ∅`) | `∂ᵢM = 0`, vanishes | gradient **exact**; `c_glob` = physical strut width |
| **varied** (`A ≠ ∅`) | dropped | controlled first-order approximation; valid where `M` varies slowly over a cell (keep `σ_m` > seed spacing). Error largest at source boundaries. |

Graded thickness enters only through the additive shift `−c_loc + c_glob`, so
the surface always sits at `c_glob` regardless of how local width varies.

### 5.5 Post-processing (order matters)

1. **Taubin smoothing** of the field — `λ|μ` band-pass (`0 < λ < −μ`), 6-neighbour
   Laplacian with per-axis weights `∝ 1/h_a²` (isotropic in world space for
   non-cubic voxels; reduces to the plain `/6` average when cubic).
   Skips voxels `>= air_skip_level()`.
2. **Boundary clamp** — `Φ ← max(Φ, Ψ_Ω + c_glob)` (Boolean intersection with the
   container). Runs *after* smoothing so the wall stays sharp.
3. **Island removal** — BFS on the solid phase, keep the largest component.
4. **Seal** — force the outer voxel shell to air so MC caps truncated features.
5. **Porosity** — counted here (§7.1).

### 5.6 Iso-surface extraction

Lewiner MC33 — topologically consistent for ambiguous cases. Buffers are sized
from two exact counting pre-passes (mirroring the real crossing tests) rather
than the theoretical worst case, which at billion-voxel grids would exhaust RAM.

---

## 6. Guards

- `seeds.size() < 3` → error, return false.
- `isoLevel <= 0` → error, return false (a non-positive iso-level produces a
  field that never crosses the surface, i.e. an empty mesh — `v ≥ 0` always).
- `compute_scalar_field` / `marching_cubes` return `bool`; every caller chains on
  success so a failed stage skips the rest instead of running on garbage.

---

## 7. Metrics

| metric | substrate | resolution |
|---|---|---|
| Volume, Surface, BS/BV, BS/TV, SMI | triangle mesh | generation grid |
| Porosity (primary) | scalar field | generation grid |
| Porosity (legacy) | mesh vs. analytic container volume | — |
| Tb.Th, Tb.Sp | resampled binary image | `voxelSize` |
| Tb.N, DA | resampled binary image | `voxelSize` |
| Conn.D (voxel) | resampled binary image | `voxelSize` |
| Conn.D (mesh, legacy) | triangle mesh | generation grid |

### 7.1 Porosity

Counted on the **generation grid**, numerator and denominator sharing one
discretization so the bias cancels:

```
P = 1 − |{x : Ψ_Ω(x) ≤ 0 ∧ Φ(x) ≤ c_glob}| / |{x : Ψ_Ω(x) ≤ 0}|
```

The **legacy** `P = 1 − V_mesh/V_Ω` is retained as `porosityMesh` but
systematically **over**estimates: the sealed, sub-voxel MC surface can never
reach the container wall. `BV/TV = 1 − P`.

### 7.1b Surface normalizations — BS/BV vs BS/TV ⚠️

Two **different** standard quantities (Parfitt/ASBMR). They differ by a factor of
`BV/TV`, i.e. roughly an order of magnitude. Both are computed and reported:

```
BS/BV = surfaceArea / volume         (surface per BONE volume)
BS/TV = surfaceArea / domainVolume   (surface per TOTAL sample volume)
BS/TV = (BS/BV) * (BV/TV)
```

`BS/BV` is tied to thickness and is a useful cross-check:
`BS/BV ≈ 2/Tb.Th` for ideal plates, `≈ 4/Tb.Th` for ideal rods.

> **Do not compare across normalizations.** When reading a reference value,
> sanity-check it against the source's own Tb.Th and BV/TV before using it as a
> target — see §9 for a case where a published "BS/TV" value is provably not
> BS/TV.

### 7.2 Local thickness / separation (Hildebrand–Rüegsegger)

Separable exact EDT (1-D chamfer along x → square → Felzenszwalb parabola along
y and z) → distance-ridge (26-neighbour sphere-inclusion test) → paint each ridge
sphere, keep the max diameter per voxel. Separation is the same routine on the
**inverted** image; out-of-container voxels are forced to background so pores are
bounded by the wall.

The mean is taken over the **foreground only** (`field != 0`). Gating on
`thicknessMap > 0` instead would leak the one-voxel shell of background that
ridge spheres over-paint at their boundary.

**Known biases (quantified by the slab phantom, §8):**
- **Measurement bias ≈ +2 voxels on the diameter**, vanishing linearly:
  `measured ≈ T + 2·voxelSize`. At `v = 0.05` that is **+82%** on a 0.123 mm
  feature; at `v = 0.01`, +13%. **Match `voxelSize` to the reference µCT.**
- **Junction bleeding**: the metric assigns each voxel the diameter of the
  largest sphere *containing* it, so thick Plateau-border junctions inflate the
  thickness of adjacent thin walls. This makes the mean Tb.Th sit *above* the
  nominal `c_glob` by real geometry — not an error. It is also why smin fillets
  (§5.3) must be avoided or kept tiny.

### 7.3 Tb.N and DA (MIL)

Fibonacci-sphere directions × a stratified grid of parallel rays, walked over the
**resampled image** at `voxelSize`.

- **Container-masked (non-box correctness).** Both MIL walks restrict the
  accumulated intercept length **and** the interface-crossing count to voxels
  whose centre is **inside the container** (the same `is_inside` mask used by
  Tb.Th/Tb.Sp/porosity). Without it, the resampled AABB of a non-box container
  (e.g. a **cylinder**) carries empty, structure-free corners *outside* the wall:
  their ray length would dilute `P_L` and bias Tb.N **low** (≈ −21% for a cylinder
  inscribed in its box, the `4/π` area ratio) and — because the corner traversal
  is *direction-dependent* (axial rays skirt them, radial rays cross them) —
  distort the DA ellipsoid. A **box** fills its AABB, so its mask is all-inside and
  Tb.N/DA are unchanged. The container wall itself (interior solid → exterior air)
  is not counted as a crossing.
- **Tb.N** = `P_L / 2 / voxelSize`. Each trabecula crossed by a test line gives
  **two** interface crossings (entry + exit), hence the `/2`. *(An earlier version
  dropped it and read 2× high.)*
- **DA** = ratio of the MIL ellipsoid's radii (fitted quadric → eigenvalues).
  Being a **ratio** it is scale-invariant in *units*, but ⚠️ **not in sampling**:
  coarser voxels blur the ellipsoid toward isotropy and pull DA toward 1
  (measured: 1.772 → 1.636 from v = 0.01 → 0.05, see §9). Always quote DA at a
  stated voxel size and calibrate the stretch at the µCT resolution.

### 7.4 Connectivity density

**Voxel Euler characteristic (preferred).** χ = V − E + F − C of the solid's
cubical complex, summed via a 256-entry LUT over 2×2×2 lattice windows (each
cell attributed to its max-corner window: `Δχ = vertex − 3 edges + 3 faces −
cube`). Out-of-range = background (implicit padding closes edge-truncated
features). `Conn = 1 − χ`; `Conn.D = Conn / TV`. One O(N) parallel pass.
Currently **6-connectivity** only (26 is a future LUT).

**Mesh surface-genus (legacy).** `genus = 1 − χ_surface/2`, `Conn.D = genus/V_Ω`.
Valid **only** for a single, cavity-free closed surface: the correct relation is
`genus = C − χ/2` for `C` components, and every enclosed void is its own surface
component. Closed-cell foams therefore break it.

Verified: both agree on cavity-free structures (open lattice: 5.09 vs 5.11;
near-solid block: ≈0), and diverge exactly where cavities appear (τ=0.3: mesh
2.95 vs voxel 2.06).

> `Conn = 1 − χ = β₁ − β₂` when cavities exist, so cavity-heavy structures can
> read low/negative. This is inherent to Euler-based connectivity (BoneJ too);
> the win is that the voxel method computes the *correct* solid χ.

### 7.5 SMI (Structure Model Index)

```
SMI = 6·(BS′·BV)/BS²        BS′ = dBS/dr
```

`BS′` by forward difference: displace every vertex along its (normalized) vertex
normal by a small `δ` (default `0.05·h`) and re-sum triangle areas. Dilating a
**plate** barely grows its area (`SMI ≈ 0`); dilating a **rod** grows its lateral
area (`SMI ≈ 3`); a sphere gives ≈4.

**Reported, not used for calibration** — see §9.

---

## 8. Validation & calibration tools

### 8.1 Analytic phantom suite (`--phantom T v`, `--phantom-suite f v`)

Shapes with **known ground truth**, built directly (bypassing the generator) via
`build_phantom_field`, to verify the estimators independently of the generator.
All are bounded in three axes (an infinite slab is degenerate for the separable
EDT — two axes have no background and the squared distance overflows).

| shape | Tb.Th | SMI | Conn (β₁) |
|---|---|---|---|
| slab (plate) | = feature | ~0 (ideal) | 0 |
| cylinder (rod) | = feature | 3 | 0 |
| sphere | = feature | 4 | 0 |
| torus (rod tube) | = feature | ~3 | 1 |

**Tb.Th convergence** (single slab, `--phantom`):

| voxel | true T | measured | bias |
|---|---|---|---|
| 0.05 | 0.123 | 0.223 | +2.0 vox |
| 0.02 | 0.123 | 0.163 | +2.0 vox |
| 0.01 | 0.123 | 0.139 | +1.6 vox |

→ the thickness estimator is **asymptotically unbiased**; bias `≈ 2·voxelSize`.

**Suite** (feature 0.3 mm, `--phantom-suite`):

| shape | Tb.Th (true 0.3) | SMI | SMI true | Conn | Conn true |
|---|---|---|---|---|---|
| sphere | 0.300 | 4.00 | 4 | 0 | 0 |
| torus | 0.301 | 2.99 | 3 | 1 | 1 |
| cylinder | 0.293 | 3.29 | 3 | 0 | 0 |
| slab | 0.293 | 2.30 | ~0 | 0 | 0 |

- **Tb.Th**: within ~1 voxel for every shape.
- **Conn**: *exact* topological integer (torus = 1, others = 0) — validates the
  voxel-Euler connectivity implementation.
- **SMI**: exact on smooth shapes (sphere 4.00, torus 2.99); cylinder +10% from
  the flat end caps; the **slab reads 2.30, not ~0**, because SMI is sensitive to
  **sharp edges** (the box-SDF plate's 90° corners inflate under normal-offset
  dilation). The trend sphere → cylinder → slab is monotone in edge sharpness.
  This is a property of SMI itself, not the estimator — and it implies our
  faceted scaffold SMI is an upper-ish estimate. Use sphere and torus as the
  clean SMI ground-truth points.

### 8.2 Thickness calibration (`calibrate_thickness`, `--calibrate-target T v`)

Secant search over the **monotonic** `c_glob → Tb.Th` mapping, regenerating and
re-measuring each step. Because it closes the loop on the *actual* output it is
automatically agnostic to openness, seed spacing, voxel size and container shape
— no fitted model, no per-container calibration.

**The target must be measured at a fixed `voxelSize`** (the µCT's), since the
`+2v` bias means "Tb.Th" is resolution-dependent. Match measured-to-measured.

**Uniform vs varied thickness.** For uniform thickness the search variable is the
iso-level. For a *varied* (graded) thickness — a spatially interpolated iso-level
over `[startThickness, endThickness]` — the iso-level is ignored, so calibration
instead scales the whole range by one multiplicative factor (preserving the
grading ratio) until the measured **mean** Tb.Th hits the target; the GUI sets the
target to the range mean `(start+end)/2`. `calibrate_thickness` returns a `bool`
(valid calibrated scaffold) and accepts an optional `onProgress` callback so the
GUI task can drive the progress bar. It leaves the mesh built, so callers only
need `estimate_metrics` afterwards (no second `compute_scalar_field`/`marching_cubes`).

**Joint calibration** (`two_knob_calibration`, `three_knob_calibration`), built
on the §8.4 sensitivity structure — each is a *structure-exploiting* solve, not a
generic optimizer:

- **`calibrate_openness`** — 1-D secant on openness → target porosity, at the
  fixed (already-calibrated) thickness. Porosity is the voxel estimate from
  `compute_scalar_field`, so no marching cubes is needed per step. Clamped to
  [0,1]; reports "unreachable" if pinned at a bound off-target.
- **2-knob (Tb.Th, porosity)** — the Jacobian is **near-triangular**
  (∂Tb.Th/∂openness ≈ 0), so this is a back-substitution: `calibrate_thickness`
  then `calibrate_openness`, with a cheap outer check for the tiny residual
  coupling (converges in one pass). *Not* Broyden/BFGS/LM — we know the structure.
  **Junction-fillet back-off:** when `smoothJunctions` is on (§5.3b), the fillet
  fattens the nodes and raises the *minimum* reachable Tb.Th; if that floor sits
  above the target, `calibrate_thickness` converges the iso-level down while Tb.Th
  stays pinned high. The 2-knob solve detects the missed target and **halves `k`
  (`smoothK`)** — recalibrating each time, up to 6 back-offs, dropping to the plain
  linear blend once `k` underflows — until the target becomes reachable. A box /
  smoothing-off run is unaffected.
- **3-knob (Tb.Th, porosity, SMI)** — thickness first (frozen), then the
  `(openness, spread) → (porosity, SMI)` block. That block is **ill-conditioned**
  (both knobs move the outputs nearly parallel, det ≈ 0.02), so instead of
  inverting a near-singular 2×2 (a plain Newton overshoots to the domain corners —
  observed) we solve it **nested**: secant on **spread** → SMI, with
  `calibrate_openness` re-run inside each step to hold **porosity exactly**. Both
  sub-problems are monotone 1-D. Consequence (a real limitation, not a solver
  defect): **SMI is only weakly controllable** — reachable only within a narrow
  band at fixed density; targets outside settle spread at a bound (best-effort),
  with porosity still hit exactly. Verified: (Tb.Th 0.127, porosity 0.842) → hit
  to ~0.1%; SMI reached to the band edge. Its thickness step shares the 2-knob
  **junction-fillet back-off** (halve `k` until Tb.Th is reachable).

Both joint solves guard the smoothing knobs with RAII (`SmoothKnobGuard`):
`smoothJunctions`/`smoothK` are restored on any failure/cancel return and kept
only on success (the backed-off `k` is part of the calibrated result). All three
return `bool` (valid scaffold) + take `onProgress`, and leave the mesh built. GUI checkboxes map to them with the constraint **SMI ⇒ porosity** (SMI
cannot be targeted alone — it needs the openness+spread pair): {thickness},
{thickness+porosity}, {thickness+porosity+SMI}. Getters `get_iso_level()`,
`get_openness()`, `get_spread()` expose the calibrated "background" knob values so
the GUI can show target vs calibrated-knob vs achieved.

### 8.3 Replicates and reproducibility

**The seed placement is stochastic.** A Poisson radius fixes the *minimum
spacing*; it does not fix *which* of the many valid arrangements at that spacing
you get. Bridson's dart-throwing draws random candidates, so every run yields a
different seed cloud — same statistics, different geometry. One generation is a
single draw from a **distribution** of statistically-equivalent scaffolds.

**Measured impact is small.** Five realizations of the Sample 1 recipe
(seeds 1..5, measured at 0.014):

| metric | mean ± SD | CoV |
|---|---|---|
| Tb.Th | 0.1254 ± 0.0008 | 0.7% |
| Tb.Sp | 0.7207 ± 0.0105 | 1.5% |
| Tb.N | 1.052 ± 0.020 | 1.9% |
| DA | 1.916 ± 0.026 | 1.3% |
| SMI | 1.456 ± 0.055 | 3.8% |
| BS/BV | 26.60 ± 0.09 | 0.3% |
| Conn.D | 4.384 ± 0.311 | 7.1% |
| porosity | 0.8869 ± 0.0020 | 0.2% |

The generator is **statistically very reproducible**: each metric averages over
~10⁶ solid voxels (or 40M MIL rays), so per-seed randomness washes out by the law
of large numbers. Conn.D is the noisiest (7%) because it counts a small integer
number of loops rather than averaging a field.

> ⚠️ **These SDs are not comparable to the literature's SDs.** Ours measures the
> *numerical* spread of the generator at **fixed** parameters. The literature's
> (e.g. DA 1.74 ± 0.15) measures **biological variability across specimens**. Ours
> being ~6x tighter does not mean we reproduce their spread — it means the
> generator is precise. To mimic the biological spread you must **sample the input
> parameters** across the reported ranges, not merely re-seed.

> ⚠️ Because the realization-to-realization spread is this small (DA ±0.026), a
> difference of ~0.2 in DA between two runs is **systematic, not stochastic** —
> look for a settings or code-version difference, not luck.

> The metrics themselves are **deterministic**: `estimate_anisotropy` and
> `estimate_trabecular_number` seed their ray jitter with
> `std::mt19937 rng(1337 + omp_get_thread_num())`, a fixed seed. All run-to-run
> variation comes from the seed generator, not the measurement.

**`InterfaceSeedGenerator::rngSeed`** controls this:
- non-zero (**default `1`**) — reproducible: regenerating at identical parameters
  gives an identical cloud (and identical metrics, incl. DA); distinct values give
  independent realizations.
- `0` — opt-in non-deterministic: a fresh `std::random_device` draw per run.

The default is a **fixed** seed so the GUI (which does not expose it) maps
identical inputs to identical output — a deterministic modelling tool. A metric
difference between two GUI regenerations at the same parameters is therefore a
real settings/version difference, **not** seed luck (the seed is now constant).
Realization spread is obtained deliberately via the profiler (`--replicates` /
`--seed`), which sets the seed explicitly.

Each `run()` now seeds **one** RNG stream, shared by the root fallback and the
candidate sampling, so a given `rngSeed` reproduces the entire cloud exactly.

**Protocol:** report **mean ± SD over N ≥ 5 realizations** (seeds `1..N`), which
is the direct analogue of the literature's mean ± SD across *specimens*. Never
report a single run. `--replicates N v` does this and prints a per-run CSV plus
the summary; `--seed N` reproduces one specific realization.

### 8.4 Parameter → metric control map (sensitivity study)

Measured with `--sweep-input` / `--sweep-grid` around the Sample-1 baseline
(box 4 mm, radius 0.684, thickness 0.127, openness 0.6, stretch (1.1,2.0,1.1),
spread 0.5), fixed seed, measurement voxel = generation voxel (0.025). Values are
**elasticities** `S = d ln(output)/d ln(input)` (script `sensitivity.py`). Sign =
direction; |S| ≈ 1 means proportional. Seed-noise floor (5 replicates): Tb.Th
±0.3%, porosity ±0.4%, DA ±1.2%, Tb.Sp ±2.3%, Tb.N ±3.0%, SMI ±4.4%, Conn.D ±6.1%.

| output ↓ / input → | thickness | openness | stretch | spread | radius |
|---|---|---|---|---|---|
| **porosity** | −0.37 | +0.15 | −0.04 | −0.02 | +0.28 |
| **Tb.Th**    | **+0.79** | **−0.02** | −0.06 | −0.01 | −0.02 |
| **Tb.Sp**    | −0.37 | +0.18 | −0.32 | −0.03 | **+1.35** |
| **Tb.N**     | +0.54 | −0.43 | +0.23 | +0.07 | **−1.71** |
| **DA**       | +0.03 | −0.06 | **+0.78** | −0.01 | +0.05 |
| **Conn.D**   | −0.74 | +1.00 | −0.09 | −0.06 | **−2.84** |
| **SMI**      | −2.90 | +2.81 | −0.53 | −0.19\* | +1.14 |
| **BS/BV**    | −1.11 | +0.26 | 0.00 | −0.04 | +0.15 |

Per-knob reading:

- **radius (Poisson spacing) → Tb.Sp** (+1.35, the primary spacing control). Being a
  *structural* input (it re-seeds the cloud), it also strongly drives everything
  tied to cell count: Tb.N (−1.71), Conn.D (−2.84), BS/TV (−1.70), SMI (+1.14).
  Crucially it leaves **Tb.Th (−0.02) and DA (+0.05) untouched**, so spacing can be
  fixed first.
- **stretch → DA** (+0.78, the clean anisotropy control). Negligible on porosity
  (−0.04) and Tb.Th (−0.06), so anisotropy can be set second without disturbing
  the density/thickness solve.
- **thickness (iso-level) → Tb.Th** (+0.79, primary). Also lowers porosity (−0.37)
  and BS/BV (−1.11, matching the exact `BS/BV ∝ 1/Tb.Th` law — a pipeline sanity
  check). So thickness is one of the two density knobs.
- **openness → porosity/BV·TV** and the mean rod↔plate axis: strong on SMI (+2.81),
  Conn.D (+1.00), Tb.N (−0.43). **Its elasticity on Tb.Th is ≈0 (−0.02)** — the key
  result: openness barely perturbs Tb.Th, so (thickness, openness) → (Tb.Th,
  porosity) is a **near-triangular** system (Jacobian below).
- **spread → SMI, at matched density.** The *raw* elasticities are all ≤ 0.19
  (\*), which made spread look inert — but that sweep held openness fixed, so
  spread's BV/TV leak masked its true effect. At **fixed porosity** (openness
  re-calibrated to compensate), spread moves **SMI +49%** (0.50→0.75 over spread
  0→1) with Tb.Th flat. Proof it is a real degree of freedom: at porosity ≈ 0.847,
  spread 0 → SMI 0.55 but spread 1 → SMI 0.82 — the *same* porosity at *different*
  SMI, which openness alone cannot reach. So **openness and spread jointly span the
  2-D (porosity, SMI) plane**: openness traces the biological perforation curve
  (SMI and BV/TV locked together), and spread moves *off* that curve.
- **Rod/plate ↔ BV/TV are coupled by construction** (a rod encloses less material
  than a plate: `wallVal ≤ strutVal`), which is physically faithful — you cannot
  change mean rod/plate at fixed density with *either* knob alone. Spread is the
  tool that breaks the coupling when you need to.

Jacobian for the joint calibrator (`sensitivity.py`, thickness×openness grid):
```
 ∂(TbTh)/∂thickness = +0.810   ∂(TbTh)/∂openness = +0.006  ≈ 0   (near-triangular)
 ∂(por )/∂thickness = −1.769   ∂(por )/∂openness = +0.220
 det = 0.189   cond = 20.3
```

**Replication workflow (5 knobs, applied in decoupling order):**
1. **Tb.Sp** → Poisson radius (seed generation).
2. **DA** → stretch (dominant-axis anisotropy).
3. Target **Tb.Th + porosity** → solve (thickness, openness). Near-triangular:
   thickness sets Tb.Th (openness ≈ irrelevant), openness then sets porosity.
4. *(3b)* If **SMI** is also a target → add **spread** as the third knob:
   solve (thickness, openness, spread) → (Tb.Th, porosity, SMI). For sources that
   do not report SMI, spread is held at a default and SMI reported as an outcome.
5. **All remaining metrics** (Tb.N, Conn.D, BS/BV, BS/TV, tortuosity) are then
   measured, not tuned.

### 8.5 Other profiler modes

- `--calibrate t0 t1 n v` — iso-level sweep, prints the `c_glob → Tb.Th` curve.
- `--calibrate-target T v` — closed-loop secant calibration to a target Tb.Th.
- `--metrics v` — one-shot: all metrics at a single measurement voxel size.
- `--sweep-voxel a b s` — resolution convergence (§12.2b).
- `--phantom T v`, `--phantom-suite f v` — analytic ground truth (§8.1).
- `--replicates N v` — N realizations, per-run CSV + mean ± SD (§8.3).
- `--seed N` — reproduce one specific realization.
- `--sweep-radius-std rMin rMax a b n` — stochastic-radius Poisson (§3): sweep the
  per-seed radius `σ` over `[a,b]` in `n` steps and print `radiusStd, seedCount,
  TbSp_mean, TbSp_intraStd`, i.e. the input-`σ` → output-Tb.Sp transfer (mean and
  intra-sample SD). Averages `--replicates` realizations (default 3).
- `--bonej-suite N v` — emit a cohort of N scaffolds for the BoneJ cross-tool
  comparison (§12.2c): each as `.stl` + `.nrrd`, plus `generation_parameters.csv`
  and `scaffoldbuilder_metrics.csv`. Requires `--box` (box-only experiment).
  `--bonej-out dir`, `--bonej-thickness-range a b`, `--bonej-spacing-range a b`.
  The `.nrrd` is exported over the **full container box** (not the mesh AABB) so
  TV matches BoneJ's — otherwise porosity/BS-TV/BV-TV would use different total
  volumes and not be comparable. Measurement voxel `v` = the source µCT size
  (from §12.2b). Warning: the `.stl` files are large (~2M triangles → ~100 MB
  each at a 0.025 generation grid); the `.nrrd` is what BoneJ consumes.

### Configuration: CLI vs. JSON

`nlohmann/json` 3.11.3 is already vendored (`include/json.hpp`, used by
`main.cpp`), so a config-file front end is nearly free. Recommendation: **keep
both, layered** — the CLI for one-off/interactive runs, and a JSON config for the
paper's batch experiments (11 samples x replicates x sweeps), where a committed
config is a far better reproducibility artifact than a shell-history line. Add
the JSON path when the multi-sample batch (§12.4a) is built, not before.

---

## 9. Empirical findings

**Spread is not an independent rod/plate knob.** (4 mm box, thickness 0.13.)

At fixed `threshold = 0.5`, raising `spread` appears to change morphology:

| spread | Tb.Th | porosity | SMI | Conn.D |
|---|---|---|---|---|
| 0.0 | 0.1357 | 0.879 | 1.16 | 3.95 |
| 0.8 | 0.1360 | 0.851 | 0.68 | 3.02 |

But at **matched porosity** the SMI difference vanishes:

| spread | threshold | porosity | SMI |
|---|---|---|---|
| 0.0 | 0.50 | 0.8784 | 1.130 |
| 0.8 | 0.60 | 0.8759 | 1.151 |

→ `spread`'s only measurable effect is on **BV/TV** (via Jensen: material is
nonlinear in τ, so spreading τ about a fixed mean nets more solid). SMI is a
*global average* and is blind to the *heterogeneity* that `spread` adds.

### Three different "SD"s — do not conflate ⚠️

| SD | measures | typical (Sample 1) |
|---|---|---|
| `Tb.Th_SD` (metric output) | spatial spread of wall thickness *within one scaffold* | ±0.021 |
| replicate SD (`--replicates`) | generator reproducibility *across seeds*, fixed inputs | ±0.0008 |
| literature SD (e.g. Ulrich) | biological variability *across specimens* | ±0.017 |

Only the last two are loosely comparable. Never plot `Tb.Th_SD` as an error bar
against a literature SD band. On a `--sweep-voxel` plot **no error bar is
appropriate at all**: each point re-measures one fixed geometry deterministically,
so the variation *is* the signal, not uncertainty.

### Resolution convergence on a real scaffold (§12.2b)

Sample 1 recipe (box 4 mm, Poisson 0.684, thickness 0.127, generation voxel
0.025, stretch 1.1/2.1/1.1, openness 0.636, spread 0.5, seed 1), generated once
and re-measured at nine image voxel sizes:

| voxel | Tb.Th | Tb.Sp | Tb.N | DA | Conn.D |
|---|---|---|---|---|---|
| 0.010 | 0.1207 | 0.750 | 1.084 | 1.772 | 3.94 |
| **0.014** (µCT) | **0.1243** | **0.761** | **1.075** | **1.762** | **3.89** |
| 0.020 | 0.1291 | 0.758 | 1.066 | 1.752 | 3.81 |
| 0.030 | 0.1388 | 0.757 | 1.052 | 1.726 | 3.83 |
| 0.040 | 0.1616 | 0.756 | 1.033 | 1.686 | 4.09 |
| 0.050 | 0.1655 | 0.756 | 1.018 | 1.636 | 4.17 |

Least-squares fits over the range (R² separates law from noise — always report it):

| metric | intercept (v→0) | R² | verdict |
|---|---|---|---|
| Tb.Th | **0.1064** | **0.9727** | real linear law |
| Tb.N | **1.099** | **0.9964** | real linear law |
| DA | **1.815** | **0.9710** | real linear law |
| Tb.Sp | (0.755 = mean; slope −0.012) | **0.0033** | **flat** — report the slope, not the intercept |
| Conn.D | — | **0.6221** | **no trend**; deliberately not fitted |

Conn.D's R² = 0.62 is the instructive one: high enough that a fitted line would
have looked plausible, low enough that it would have been fiction. It derives
from an **integer** Euler characteristic, so it changes in discrete jumps as
loops appear/vanish rather than varying smoothly.

**Why Tb.Sp is flat and Tb.Th is not.** Not because of a different measurement
path — they are the *same* algorithm on the *same* resampled image
(`estimate_local_thickness`, `separation=true` merely inverts the phase). The
difference is **feature size relative to the voxel**: pores (~0.75 mm) span
15–75 voxels over the tested range and are well resolved throughout, whereas
trabeculae (0.127 mm) span only **2.5–12 voxels**. The discretisation bias is a
**thin-feature effect**.

⚠️ The Tb.Th fit slightly **over-extrapolates**: the curve bends below the line
at the fine end (0.010–0.015), so the true asymptote is a little above 0.1064.
Either state that the linear fit is an approximation over the full range, or fit
only over the fine half where the tool actually operates.

**At the µCT resolution (0.014)** — the only fair comparison point, since the
reference measured at 14 µm:

| metric | ours @0.014 | Ulrich (mean ± SD) | inside band? |
|---|---|---|---|
| Tb.Th | 0.124 | 0.127 ± 0.017 | ✅ |
| Tb.Sp | 0.761 | 0.684 ± 0.109 | ✅ |
| DA | 1.762 | 1.74 ± 0.15 | ✅ |
| Tb.N | 1.075 | 1.45 ± 0.20 | ❌ (the structural gap) |

Three of four land inside the reference's own SD band; Tb.N is the known
Voronoi-vs-rod-and-plate topology gap. Note the distinction worth stating in the
paper: the *resolution-free* DA of this geometry is 1.81, but we measure 1.762 at
14 µm and the reference measured 1.74 at 14 µm — **we match the reference at the
reference's resolution**, which is the correct comparison.

- **Controls held exactly**: SMI 1.21334, BS/BV 25.8523, BS/TV 3.00422, porosity
  0.882441 — *bit-identical* at every step, confirming the sweep isolates
  measurement resolution and no state leaks between steps.
- **Tb.Th ≈ T + k·v** with slope ≈ 1.1 (extrapolating to T ≈ 0.110 at v→0). The
  phantom's slope was ≈2 voxels on a *plate*; a foam sits lower because junction
  bleeding partly saturates. The linear-in-`v` law holds on real scaffolds.
- **Tb.Sp is resolution-immune** (±0.6%): pores are 50–75 voxels across even at
  the coarse end. The bias is a **thin-feature** effect, not a general one.
- **Conn.D noisy but trendless** (±5%, no drift) — it counts a small integer of
  loops rather than averaging a field.

⚠️ **DA drifts with resolution: 1.772 → 1.636 (−8%) as voxels coarsen.** DA is
scale-invariant in *units* (it is a ratio) but **not** in *sampling*: coarser
voxels blur the MIL ellipsoid toward isotropy, pulling DA toward 1. So DA must be
**quoted at a stated voxel size and calibrated at the µCT resolution** — tuning
stretch at 0.025 leaves it ~1.5% low at 0.014. At the µCT-matched 0.014 this
recipe gives DA ≈ 1.765, matching Ulrich's 1.74 ± 0.15 (stretch 2.1; the earlier
2.6 was over-corrected).

**Consequences:**
- Do **not** use `spread` as an SMI knob in calibration — it would fight openness
  for BV/TV. Tune BV/TV with `threshold` alone.
- `spread` remains a **structural-realism** option: rod-and-plate coexistence is
  geometrically real and visible even though SMI cannot score it.
- SMI is **reported**, not targeted. Independent rod:plate control would need a
  *heterogeneity* descriptor (per-trabecula/ITS-style classification), not global SMI.

**`spread` is orthogonal to thickness** — Tb.Th moves <0.2% for spread 0→0.8, so
a calibrated `c_glob` survives spread changes.

**Tb.N gap vs. literature is structural, not a bug.** Even the reference's own
numbers imply it: `1/(Tb.Th+Tb.Sp) ≈ 1.23` while they report 1.45 — real bone is
a mixed rod-and-plate network with extra branch nodes, so it crosses more
trabeculae per mm than a Voronoi foam of matched cell size.

More precisely, for Ulrich Sample 1:
- **Reference**: Tb.Th+Tb.Sp = 0.811 → plate model predicts **1.23**, they report
  **1.45** (+18% *above* the model: extra branch nodes).
- **Ours**: Tb.Th+Tb.Sp = 0.830 → model predicts **1.20**, we measure **1.03**
  (−14% *below*).

The Voronoi foam has *fewer* trabeculae than its own spacing implies; real bone
has *more*. No parameter closes this — it is a topology difference.

### Ulrich (1999) Table 1 "BS/TV" is actually BS/BV ⚠️

Ulrich et al. (1999), *Bone* 25(1):55–60, Table 1 reports (units as printed):

| site | BV/TV (%) | "BS/TV" (/mm) | Tb.Th (mm) | Tb.Sp (mm) | Tb.N (/mm) | DA |
|---|---|---|---|---|---|---|
| Calcaneus | 11.65 | 21.82 | 0.127 | 0.684 | 1.45 | 1.74 |
| Femoral head | 20.67 | 15.25 | 0.172 | 0.706 | 1.42 | 1.73 |
| Iliac crest | 15.22 | 17.63 | 0.150 | 0.754 | 1.39 | 1.50 |
| Lumbar spine | 8.15 | 23.37 | 0.123 | 0.800 | 1.26 | 1.43 |

The paper's *methods* define BS/TV correctly (BS normalized by total volume), and
the column says `/mm`. **But the tabulated values are BS/BV.** Proof — the
product `Tb.Th x "BS/TV"`:

| site | Tb.Th | "BS/TV" | product | BV/TV |
|---|---|---|---|---|
| Femoral head | 0.172 | 15.25 | **2.62** | 20.67% |
| Iliac crest | 0.150 | 17.63 | **2.64** | 15.22% |
| Calcaneus | 0.127 | 21.82 | **2.77** | 11.65% |
| Lumbar spine | 0.123 | 23.37 | **2.87** | 8.15% |

The product is essentially **constant (2.62–2.87)** and lies inside `[2, 4]` —
the plate→rod range of `BS/BV · Tb.Th` (`BS/BV ≈ 2/Tb.Th` plates, `4/Tb.Th` rods).

Were these genuinely BS/TV, the product would equal `k·(BV/TV)` and would have to
swing **2.5x** across these sites (0.22 at lumbar's 8.15% up to 0.56 at femoral
head's 20.67%). It varies by 10%. Therefore the column is **BS/BV**.

A second, independent check on the calcaneus row: as BS/TV it would imply
`BS/BV = 21.82/0.1165 = 187 /mm` → plate-model `Tb.Th = 2/187 = 0.011 mm`, a 12x
contradiction of their own reported 0.127.

**Consequence:** compare our **BS/BV** against their tabulated value.

| | reference (calcaneus) | ours | Δ |
|---|---|---|---|
| BS/BV | 21.82 ± 3.09 | 26.80 | +23% (1.6 SD) |

Applying the same product test to our scaffold gives `0.125 × 26.8 = 3.35` vs
their `2.77` — closer to the rod limit (4), **independently confirming our
structure is more rod-like**, in agreement with SMI (1.5–1.7) and the Tb.N
deficit. Three metrics, one root cause.

**Lesson:** always sanity-check a reference surface value against the source's own
Tb.Th and BV/TV before adopting it as a target. Reporting **both** BS/BV and
BS/TV is what made this detectable.

### Sample 1 replication (Ulrich, healthy calcaneus)

Inputs: box 4x4x4, thickness 0.127, generation voxel 0.025, Poisson 0.684,
stretch (1.1, 2.6, 1.1), openness 0.77, spread 0.5; measured at voxel 0.014.

| metric | literature | ours | Δ |
|---|---|---|---|
| Tb.Th | 0.127 ± 0.017 | 0.130 ± 0.019 | +2.4% ✅ |
| Tb.Sp | 0.684 ± 0.109 | 0.700 ± 0.162 | +2.3% ✅ |
| BV/TV | 11.65 ± 3.33 | 11.35 | −2.6% ✅ |
| Porosity | 88.35 | 88.65 | +0.3% ✅ |
| DA | 1.74 ± 0.15 | 1.747 | +0.4% ✅ |
| BS/BV (their "BS/TV") | 21.82 ± 3.09 | 26.80 | +23% ~ |
| Tb.N | 1.45 ± 0.2 | 1.034 | −29% ❌ |
| Conn.D | not reported | 4.55 | — |
| SMI | not reported | 1.73 | — |

Five of seven within ~3%, all inside the reference's own SDs. Tb.N is the one
clear miss and is the structural limitation above.

**Cross-check:** BS/BV = 26.8 lies between the plate (15.4) and rod (30.8) limits
for Tb.Th = 0.130, leaning rod — consistent with the measured SMI of 1.73. Two
independent metrics agreeing that the structure is mid-to-rod-like.

**Stochastic spread between realizations is small** (see §8.3): DA ± 0.026,
Tb.Th ± 0.0008 over 5 realizations. Report mean ± SD anyway (it is cheap and
demonstrates reproducibility), but do **not** attribute run-to-run differences
larger than ~2% to seeding luck — those are systematic.

⚠️ **The single-run DA of 1.747 above does not reproduce.** Five replicates of
the identical recipe give **DA = 1.916 ± 0.026**, i.e. the 1.747 sits ~6 SD away
and cannot be a stochastic draw. It predates the refactor that moved anisotropy
onto the resampled image at the target voxel size (§7.3), so it was measured on
the generation grid instead. **Re-run all samples on the current build**; the
apparent bullseye on DA (1.747 vs 1.74) is likely an artifact of the old code
path, and at stretch (1.1, 2.6, 1.1) the current DA of 1.92 is ~1.2 SD *above*
the reference's 1.74 ± 0.15, so the stretch needs re-tuning.

---

## 10. Known limitations / open items

- OpenMP `collapse(3)` is ignored under MSVC's OpenMP 2.0 (warning C4849) —
  correctness unaffected (loops still parallelize over the outer index); only the
  coarse narrow-band loops lose some parallelism.
- 26-connectivity Euler LUT not implemented (6 only).
- `estimate_trabecular_number` formula 1 ("Derived Proxy") **silently falls back
  to MIL** if thickness has not been re-estimated since the last Tb.N run.
- Varied-anisotropy kNN uses a Euclidean candidate set (`3κ²` heuristic).
- SMI: global average only; cannot see morphological heterogeneity.
- The GUI does not expose `rngSeed`. It now defaults to a **fixed** value (`1`),
  so GUI scaffolds are reproducible (identical inputs → identical output); the
  cost is that the GUI cannot draw *different* realizations at fixed parameters —
  use the profiler (`--replicates` / `--seed`) for that.

---

## 11. Changelog

- **2026-07-26c** — **Reproducible-by-default seeding.** `InterfaceSeedGenerator::
  rngSeed` default flipped `0 → 1` (fixed). The GUI does not expose the seed, so it
  was re-drawing `random_device` every regeneration → metrics (notably DA, ±0.026)
  wobbled at identical parameters. With a fixed default the GUI maps identical
  inputs to identical output. Verified the DA *measurement* is already
  deterministic (two identical profiler runs → bit-identical DA per seed); the
  wobble was purely the seed cloud. `rngSeed = 0` remains an opt-in for random
  draws; the profiler still sets seeds explicitly for `--replicates`/`--seed`.
- **2026-07-26b** — **Stochastic-radius Poisson** (§3). New varied-Poisson mode:
  each seed draws its minimum spacing from a truncated normal `N(μ, σ)` on
  `[rMin, rMax]` (no distance field), from the seeded RNG stream so the cloud stays
  reproducible. Injects intra-sample spacing heterogeneity for realism. Verified
  and **quantified the non-triviality**: input `σ` ≠ output Tb.Sp SD (packing
  regularises variance; raising `σ` mostly shifts the *mean* denser). Added the
  profiler map `--sweep-radius-std` (§8.5) to read the transfer, and GUI radio
  "Random (Normal)" with μ/σ inputs on both varied-Poisson panels. Fixed the
  inverted `@param k` doc on `smin_gradient`. *Not persisted by
  `export_parameters` yet (μ/σ/flag) — follow-up.*
- **2026-07-26** — **Container-masked MIL + junction fillet documented.**
  (1) Both MIL walks (Tb.N and DA) now restrict intercept length and crossing
  count to voxels **inside the container** (§7.3), fixing a cylinder bias
  (Tb.N ≈ −21%, DA ellipsoid distortion) from the empty AABB corners; box results
  unchanged. (2) The optional **junction fillet** (`smoothJunctions`/`smoothK`,
  Quilez `smin` in `compute_scalar_field`) is now described in §5.3b — a fillet of
  radius `k` between the linear-blend field and a shifted strut field, rounding the
  rod–plate nodes; off by default. (3) `two_knob_calibration` **backs `k` off**
  (halving, ≤6×) when the fattened floor puts the Tb.Th target out of reach (§8.2).
  (4) `export_parameters`/`read_parameters` extended with `Spread`,
  `BackgroundWeight`, `SmoothJunctions`, `SmoothK`, `TransitionDistance` (appended
  columns, back-compatible); the multi-source anisotropy list is still not
  persisted (belongs in the `.scaf` format). (5) GUI **Cancel** button on both
  generate popups → cooperative cancellation of the worker task.
- **2026-07-22e** — **Joint calibration** implemented (`calibrate_openness`,
  `two_knob_calibration`, `three_knob_calibration`; §8.2). 2-knob = triangular
  back-substitution (thickness→openness). 3-knob = thickness then a **nested**
  spread→SMI secant with openness re-solved inside to hold porosity exactly — the
  (openness,spread)→(porosity,SMI) block is ill-conditioned, so a direct 2×2
  Newton/LM overshoots (observed); the nested 1-D form is robust. SMI is only
  weakly controllable (narrow band) — targets outside settle at a spread bound,
  best-effort. Verified via profiler `--two-knob`/`--three-knob`: (Tb.Th 0.127,
  porosity 0.842) hit to ~0.1%. GUI mapping: SMI⇒porosity; getters added for the
  calibrated knob values. Note: my initial "no LM/Broyden needed" held for the
  well-conditioned 2-knob block but NOT the ill-conditioned 3-knob block — testing
  corrected the design (nested beats both).
- **2026-07-22d** — Sensitivity results + **spread vindicated** (§8.4). Full
  elasticity matrix measured; the (thickness, openness)→(Tb.Th, porosity) Jacobian
  is near-triangular (∂Tb.Th/∂openness ≈ 0). Key finding: at **matched porosity**
  (spread×openness grid + constant-porosity slice), spread moves **SMI +49%** with
  Tb.Th flat — spread and openness jointly span the (porosity, SMI) plane, so
  spread is a genuine rod/plate degree of freedom, not the inert knob the raw
  sweep suggested (that sweep's fixed openness let spread's BV/TV leak mask it).
  Replication workflow formalized as 5 ordered knobs.
- **2026-07-22c** — Input/output **sensitivity study** infrastructure (precursor
  to the joint thickness+porosity calibrator). Profiler modes `--sweep-input
  <name> a b n` (one-at-a-time) and `--sweep-grid n1 a b s n2 a b s` (2D coupling
  grid) over `thickness | openness | stretch | spread | radius`, measuring the
  full metric set per config at a fixed seed (measurement voxel = generation
  voxel). `doc/paper/scripts/sensitivity.py` computes normalized-sensitivity
  (elasticity) matrices + heatmap, OAT panels, and the 2×2 `(thickness,openness)
  →(TbTh,porosity)` Jacobian for the optimizer. Also fixed a **reproducibility
  bug**: `Poisson3D::run` selected the active-list index with the unseeded global
  `rand()`, so the seed cloud varied run-to-run even at a fixed `rngSeed` (seed
  count 136/141/144 → now a stable 150). Now drawn from the seeded `gen` stream —
  this also tightens the `--replicates`/reproducibility claims (§8.3).
- **2026-07-22b** — `calibrate_thickness` now: (1) supports **varied thickness**
  by scaling `[startThickness, endThickness]` by one multiplicative factor to hit
  the target **mean** Tb.Th (previously it tuned the ignored iso-level → no-op for
  graded thickness); (2) returns **`bool`** (valid calibrated scaffold); (3) takes
  an optional `onProgress` callback for the GUI progress bar. The scaffold-factory
  task now calls only `estimate_metrics` after calibration (the redundant second
  `compute_scalar_field`/`marching_cubes` removed — calibration already leaves the
  mesh built) and maps calibration progress into [0.05, 0.95]. Factory varied
  target fixed from `|start−end|·0.5` to the true mean `(start+end)·0.5`. Uniform
  path verified (target 0.127 → 0.1266). Profiler `--calibrate-target` updated for
  the new bool return.
- **2026-07-22a** — Anisotropy: **dropped DA formula mode 1** (`1 − Lmax/Lmin`),
  which is always ≤ 0 since `Lmax/Lmin ≥ 1` — not a meaningful DA. Remaining modes:
  0 = `Lmax/Lmin` (ratio ≥1, isotropic=1, matches Ulrich/ratio-convention papers),
  2 = `λmin/λmax` ([0,1], isotropic=1), 3 = `1 − λmin/λmax` ([0,1], isotropic=0,
  BoneJ). Also fixed two **inverted GUI radio-button labels** (mode 2 said
  "MaxEig/MinEig" but computes min/max; mode 3 said "MaxEig/MinEig"). Verified
  container-shape correctness of the metrics: `get_volume()` is exact per shape
  (cylinder πr²h, abstract signed-tetrahedron mesh volume) and voxel porosity
  masks by the container SDF (`containerDistField ≤ 0`), so BS/TV, BS/BV, Conn.D
  and porosity are correct for cylinder and abstract containers (the sample 10–11
  porosity gap is an openness-range tuning limit, not a bug). Documented the Tb.N
  definition gap (see §10 limitations): reference studies use different Tb.N
  conventions, and ours reads low against the direct-3D values of the high-porosity
  studies — a measurement-definition effect, not a structural mismatch.
- **2026-07-21e** — `--runtime-suite` now records **peak resident memory**
  (`peak_ram_mb`, Windows `GetProcessMemoryInfo` peak working set) per resolution
  in `runtime_vs_resolution.csv`. Headline figures at thickness 0.3 mm, 8 threads,
  5 repeats: resolution scales field N^0.96 / MC N^0.80 / total N^0.89; a µCT-grid
  scaffold (0.025 mm, 8M vox) generates in ~9.5 s using ~1 GB, the finest grid
  (0.01 mm, 125M vox) in ~126 s using ~8.5 GB; runtime is only weakly seed-count
  dependent (total n^0.29). Implementation note: `windows.h` is included with
  `NOMINMAX`/`WIN32_LEAN_AND_MEAN` and `#undef ERROR` (wingdi's `ERROR` macro
  collides with `LogPriority::ERROR`).
- **2026-07-21d** — Added `--runtime-suite` (§12.3c): benchmarks generation time
  vs grid resolution (fixed seeds) and vs seed count (fixed grid), median over
  `--runtime-repeats`, splitting scalar-field / marching-cubes / seed-gen stages
  and logging thread count. Outputs two CSVs under `--runtime-out`. Also fixed a
  seed-visualization bug: `Poisson3D::run` shadowed the member `bounds` with a
  local, so `update_model()` scaled the seed spheres from uninitialized bounds
  (huge or invisibly tiny, intermittently); now assigns the member, and `struct
  Bounds` got default `{0}` initializers as defense-in-depth.
- **2026-07-21c** — `--bonej-suite` upgrades for the redone comparison:
  `--bonej-stretch-range a b` (default 1.0–2.5) samples a per-scaffold dominant-y
  anisotropy stretch so the cohort spans a real DA range (a fixed `--stretch`
  overrides); DA is now stored in BoneJ's convention directly (`estimate_anisotropy`
  mode 3 = `1 − λmin/λmax`), removing the post-hoc square conversion. `Surface Area`
  re-enabled in `bonej_batch.groovy` (works in non-headless batch mode), and
  `bonej_comparison.py` gained an auto-detected BS comparison (resolution-dependent
  caveat: ours = MC on the continuous field, BoneJ = MC on the binary image).
  **Final agreement (n=20, redone cohort spanning DA 0.30–0.73)**: Tb.Th R²=1.00
  bias +2.9%, Tb.Sp R²=1.00 bias −2.8%, DA R²=1.00 CCC=1.00 bias +0.5%
  (nRMSE 0.9%, vs 8.9% on the earlier near-isotropic cohort — the varied-stretch
  cohort is what makes DA a real validation), Conn.D and BV/TV CCC≈1.0
  (same-algorithm implementation checks). BoneJ Volume Fraction now populates
  (BV/TV/BV columns) in non-headless batch mode; BoneJ Surface Area still emits no
  SharedTable row (BS taken from our export). `bonej_comparison.py` also made
  encoding-robust (Fiji writes Windows-1252, e.g. the `mm²` header) and label-column
  robust (SharedTable names it `\`).
- **2026-07-21b** — BoneJ agreement analysis (`doc/paper/scripts/bonej_comparison.py`,
  outputs to `doc/paper/experiments/bonej_comparison/`). Pairs our
  `scaffoldbuilder_metrics.csv` against BoneJ `bonej_metrics.csv` by scaffold id
  and reports, per metric: Pearson r/R², Spearman ρ, Lin's CCC, bias (%), RMSE
  (%), Bland–Altman LoA, and Holm-adjusted Wilcoxon signed-rank. Result over 20
  scaffolds (measurement voxel 0.025): Tb.Th R²=1.00 bias +3.4% (ours high —
  junction bleeding), Tb.Sp R²=1.00 bias −2.6%, DA CCC=0.98 bias ≈0, Conn.D and
  BV/TV CCC≈1.0 (same algorithm on same image — implementation check, not
  independent physics). **DA required formula harmonization**: BoneJ reports
  `DA = 1 − D1/D3` (eigenvalue form, [0,1]); our cohort used `estimate_anisotropy`
  mode 0 = `√(λmax/λmin)` (radii ratio ≥1). Comparing them raw is off by a square
  (gave a spurious CCC=0.36); converting our mode-0 ratio into BoneJ's definition
  (`1 − 1/DA²`) via BoneJ's exported D1/D3 columns fixes it (CCC 0.36→0.98). Lesson
  for the paper: fix ONE DA convention and state it. Reporting guidance:
  lead with bias%/Bland–Altman, treat the "significant" Wilcoxon on 2–3% biases
  as practically negligible.
- **2026-07-21a** — BoneJ automation finalized (`scripts/bonej_batch.groovy`).
  Root causes pinned down: a `.ijm` macro cannot capture BoneJ2 output (results
  live in the static `SharedTable`, whose window never appears headless); and
  `--headless` breaks Connectivity (legacy `imp.duplicate()` → `ip is null`) and
  silences Thickness/Anisotropy reporting. Fix: a Groovy script run **non-headless
  with `Interpreter.batchMode = true`**, reading `SharedTable` directly. BV/TV is
  computed in-script (`StackStatistics.mean/255`); Volume Fraction / Surface Area
  are dropped (dialogs never resolve). Comparison measurement voxel set to
  **0.025 mm** (cross-tool agreement needs matched images, not µCT resolution).
  Doc §12.2c rewritten. See also GeneratorLewiner `export_nrrd`.
- **2026-07-15e** — BoneJ comparison cohort (`--bonej-suite N v`): N scaffolds in
  a 5×5×5 box, LHS/uniform over (thickness, spacing) with an `sp ≥ 2.5·th` merge
  guard, each exported as `.stl` + `.nrrd` (`.nrrd` over the full container box
  so TV matches BoneJ) plus `generation_parameters.csv` and
  `scaffoldbuilder_metrics.csv`. Reproducible (seeded parameter sampling + seeded
  seed placement). `--sweep-extra v` added so the µCT resolution is a measured
  point in the convergence sweep.
- **2026-07-15d** — Analytic phantom suite (`build_phantom_field`:
  slab/cylinder/sphere/torus + `--phantom-suite`) verifying Tb.Th, SMI and Conn
  against known ground truth (§8.1). Resolution-convergence mode
  (`--sweep-voxel`, generate-once/re-measure) for §12.2b, with
  `--mil-directions`/`--mil-lines` to control sweep cost. Experiment design
  added as §12.
- **2026-07-15c** — Seed-placement RNG made controllable
  (`InterfaceSeedGenerator::rngSeed`, `set_rng_seed`); each `run()` now seeds a
  single stream instead of constructing `random_device` at 5 separate sites.
  Profiler gained `--seed N` (reproducible run) and `--replicates N v` (N
  realizations, per-run CSV + mean ± SD). Corrected the Ulrich BS/TV finding:
  the units really are /mm, but the product test across all four sites proves
  the column is BS/BV (the earlier 1/cm hypothesis was wrong).
- **2026-07-15b** — `BS/TV` added alongside `BS/BV` (they differ by BV/TV; a
  published reference value was found to be mislabeled — see §9). Both exported
  to CSV, shown in the GUI with tooltips. SMI added to CSV. Conn.D GUI wired
  (both algorithms write one member, so import/export stays coherent). Fixed the
  connectivity popup clearing the wrong flag, and the SMI dilation input step
  (0.1 → 0.001, `0 = auto`, clamped non-negative).
- **2026-07-15** — SMI added (reported only; shown not to track `spread`
  independently of BV/TV). Per-face fenestration (`spread`) added; dead binary/
  `mixtureSoftness` path removed; smin dropped from the blend. Voxel Euler
  Conn.D added alongside the legacy mesh method. Tb.N `/2` fixed. Anisotropy and
  Tb.N moved onto the resampled image at the target `voxelSize`. Local thickness
  mean gated on foreground. Thickness calibration + slab phantom added. Roll
  composition, scale-relative bands, porosity, and AABB/bounds fixes.

---

## 12. Experiment design (paper results structure)

Organized as a **verification → validation** hierarchy so that measurement error,
generation error, and real-bone fidelity are never conflated (the recurring
"is this a model error or a measurement error?" question). Each experiment names
the tool that produces it.

### 12.1 Capabilities showcase (qualitative)
Visual figures: the three container types (box, cylinder, abstract mesh); the
three sampling modes (random, uniform Poisson, varied Poisson); anisotropy
sources; rod↔plate morphology (openness) and fenestration (spread).

### 12.2 Verification — measurement (is the metric correct?)
Ground truth known by construction; the generator is not involved.

- **12.2a Analytic phantoms** (`--phantom-suite f v`). Recover Tb.Th, SMI, Conn
  on slab/cylinder/sphere/torus (§8.1). Establishes the estimators are unbiased
  and quantifies the `+2·voxel` thickness bias and SMI's edge sensitivity. This
  is the *strong* form of measurement validation — absolute, not cross-tool.
- **12.2b Resolution convergence** (`--sweep-voxel min max step`). **Generates
  once, then re-measures** the same geometry at each image voxel size —
  regenerating per step would confound generation randomness with measurement
  resolution. Motivates matching the measurement voxel to the µCT resolution and
  shows the dependence is understood, not hidden.
  - *Varying* (image-voxel dependent): Tb.Th, Tb.Sp, Tb.N, Conn.D.
  - *Flat* (mesh / generation-grid quantities, act as **controls**): SMI, BS/BV,
    BS/TV, porosity. A non-flat control would signal a bug.
  - Useful range **0.01–0.05 mm**. Past ~0.06 a 0.127 mm wall is under 2 voxels
    and the metrics are meaningless; below ~0.01 the cost explodes (the box is
    400³ voxels and the MIL walk dominates).
  - `--mil-directions` / `--mil-lines` trade sweep runtime against MIL sampling
    (Tb.N/DA are averages over millions of rays, so reduced counts cost little).
  - `--sweep-extra v` injects an explicit voxel size into the range. **Always add
    the µCT resolution this way** so the comparison point is a real measurement
    rather than an interpolation between sweep steps.
  - Data: `doc/paper/experiments/resolution_convergence.csv` (regenerate with the
    command in its header comment; never hand-transcribe rounded values).
- **12.2c Cross-tool agreement with BoneJ** (`--bonej-suite N v`). The *weak*
  (consistency) form: BoneJ and our pipeline are both voxel-based Hildebrand, so
  agreement proves consistency, not correctness. Complements 12.2a; does not
  replace it. A cohort of N scaffolds is sampled (LHS/uniform over thickness ×
  spacing) in a **5×5×5 mm** box, each exported as `.nrrd` (over the full box,
  TV = container) for BoneJ + `.stl` + our metrics. Measurement voxel = the µCT
  size fixed by §12.2b. Compare the **intensive** metrics (Tb.Th, Tb.Sp, Tb.N,
  DA, SMI, BS/BV) directly; the TV-normalized ones (porosity, BS/TV, BV/TV) are
  comparable only because the `.nrrd` uses the container TV.
  The measurement voxel for this comparison is **0.025 mm** (the `--bonej-suite`
  default, 200³ over a 5 mm box), *not* the finer µCT size. Cross-tool agreement
  only requires both tools to read the *same* image; it does not require µCT
  resolution (that constraint belongs to §12.4, the literature comparison). The
  coarser voxel also keeps BoneJ's local-thickness tractable (0.014 → 358³ is
  ~5× the voxels and grinds under a small heap).
  Reproduce (generate the cohort, then batch-analyze in BoneJ):
  ```
  build/scaffoldProfile --box --box-size 5 --voxel-size 0.025 --openness 0.5 \
    --spread 0 --bonej-out data/bonej_comparison --bonej-suite 20 0.025

  # BoneJ side: Fiji + BoneJ2 via a Groovy script (NOT a .ijm macro, NOT
  # headless). Writes bonej_metrics.csv into the same folder.
  Start-Process -Wait C:/Users/Kostas/Fiji/fiji-windows-x64.exe -ArgumentList `
    '--console','--run','.../scripts/bonej_batch.groovy', `
    "dir='.../data/bonej_comparison',voxel=0.025"
  ```
  BoneJ-headless caveats that dictate the above (learned empirically):
  - **A `.ijm` macro cannot capture headless output.** BoneJ2 wrappers append to
    the static `org.bonej.utilities.SharedTable` (the backing store of the
    "BoneJ results" window). Headless never displays that window, so a macro has
    no window to `selectWindow()`/save. The Groovy script reads `SharedTable`
    directly instead.
  - **`--headless` breaks two analyses outright.** Connectivity delegates to the
    legacy plugin's `imp.duplicate()` and throws `ip is null`; Thickness and
    Anisotropy run but report nothing to `SharedTable`. Running the normal GUI
    build with `ij.macro.Interpreter.batchMode = true` (set inside the script)
    suppresses windows while giving the plugins the display context they need.
  - **Volume Fraction / Surface Area dialogs never resolve** (headless or batch),
    so they are dropped. **BV/TV is computed directly in the script** as
    `StackStatistics.mean / 255` — identical to Volume Fraction (Voxel). **BS**
    is taken from ScaffoldBuilder's own export (BoneJ BS unavailable here).
  - BoneJ has no native Tb.N — derive `(BV/TV)/Tb.Th` in post-processing or
    compare our Tb.N to the literature only.
  So from BoneJ we get Tb.Th, Tb.Sp (Thickness), Conn.D (Connectivity), DA
  (Anisotropy), and BV/TV (image); BS and Tb.N come from our side.

### 12.3 Verification — generator (does it produce what you specify?)
- **12.3a Reproducibility.** `--replicates N v`: mean ± SD over N ≥ 5 seed
  realizations at fixed parameters, plus `--seed` determinism. Frame as
  *numerical* precision, NOT biological spread (§8.3 warning).
- **12.3b Controllability + calibration.** Each input knob monotonically drives
  its target metric (openness→BV/TV via `--calibrate`; spacing→Tb.Sp;
  stretch→DA; thickness→Tb.Th). Closed-loop `calibrate_thickness` hits a target
  Tb.Th. This is the "input–output consistency" experiment, reframed: report it
  as calibrated target-matching (`c_glob` is not Tb.Th, it is calibrated to it).
- **12.3c Performance / runtime** (`--runtime-suite`). Two sweeps establishing
  practicality, each timed `--runtime-repeats` times with the MEDIAN reported
  (robust to OS scheduling noise): (A) generation time vs grid resolution at a
  fixed seed cloud (voxel 0.05→0.01 mm), splitting `compute_scalar_field` vs
  `marching_cubes`; (B) generation time vs seed count at a fixed grid,
  additionally timing seed generation and recording surface area + porosity.
  The Sweep-B Poisson radii are set to multiples {2,2.5,3,4,5,7}× the wall
  thickness (isoLevel), i.e. **spacing ≥ 2× thickness**, so every structure
  stays an OPEN lattice. This matters: below ~2× thickness the walls merge into
  a near-solid, the iso-surface area collapses, and marching-cubes time becomes
  *non-monotonic* in seed count (it tracks surface area, not seed count — an
  empirically confirmed structural artifact, not a scaling law). Reports OpenMP
  thread count. Writes `runtime_vs_resolution.csv` and `runtime_vs_seeds.csv` to
  `--runtime-out` (default `doc/paper/experiments/runtime`). Run on a box:
  `scaffoldProfile --box --box-size 5 --poisson --thickness 0.3 --runtime-suite`.

### 12.4 Validation — real bone (does it match reality?)
- **12.4a Multi-site replication.** 11 specimens from three µCT studies
  (Ulrich 1999: calcaneus, femoral head, iliac crest, lumbar spine; + two more
  studies). Feed the reported Tb.Th/Tb.Sp; report our full metric set as
  mean ± SD over replicates against the literature mean ± SD. Spans BV/TV
  ~8–21%, DA ~1.4–1.9.
- **12.4b Abstract container + ROI.** Generate inside a femoral-head mesh
  (a validated recipe from 12.4a); extract an ROI; verify the ROI metrics match
  the full-scaffold metrics on a homogeneous region (ROI-consistency check).

### 12.5 Discussion / limitations
Tb.N structural gap (Voronoi foam vs. rod-and-plate bone, §9); measurement bias
and its resolution dependence; SMI blindness to morphological heterogeneity;
BS/BV vs BS/TV nomenclature (§9). Report Conn.D and SMI as descriptors where no
reference value exists.

### Minimal core (if space-constrained, e.g. IEEE brief)
12.2a phantoms · 12.3a reproducibility · 12.3b calibration · 12.4a multi-site.
This answers: *the measurement is correct, the generator is controllable and
reproducible, and it spans the physiological range.*

### ⚠️ Before building any table
Re-run every sample on the current build. Tb.N (`/2`), Tb.Th (foreground gate),
and DA (resampled-image refactor) have all changed since earlier numbers.
