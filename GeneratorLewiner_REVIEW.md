# GeneratorLewiner / Anisotropy review — 2026-07-09

Scope: `include/ScaffoldGenerator/GeneratorLewiner.cpp`, `GeneratorLewiner.h`,
`include/ScaffoldGenerator/Anisotropy.h`, `include/Utils/Utils.cpp`
(`rotation_from_direction_roll`), `include/SeedGenerator/Container.h`.

---

## 1. Anisotropy: the M tensor and the covariance formulation

### 1.1 What the code actually implements (this part is consistent — but misnamed)

`AnisotropySource::M` is **not** the Riemannian metric tensor. It is the
covariance-like *structure tensor*

    C = Rᵀ S² R        (Anisotropy.h:138, :178)

`blend_metric()` linearly blends these C tensors with Gaussian weights and
returns `Mblend.inverse()` (Anisotropy.h:254). Since R is orthogonal,
`(Rᵀ S² R)⁻¹ = Rᵀ S⁻² R`, so the *effective* metric used in
`aniso_distance_sq` is

    M(x) = C(x)⁻¹  ≈  Rᵀ S⁻² R      (single dominant source)

This is exactly consistent with the uniform-anisotropy path in
`compute_scalar_field` (GeneratorLewiner.cpp:346-355), which warps seeds by
`S⁻¹R(x − c)`, i.e. squared distance `δᵀ Rᵀ S⁻² R δ`. Good — the two paths
agree. But rename things for the paper and the code: the stored tensor is a
covariance C; the metric is M = C⁻¹. The comment `// M = R S^2 R^T` at
Anisotropy.h:177 also contradicts the code (`Rᵀ S² R`).

### 1.2 BUG — roll angle composed on the wrong side (`Utils.cpp:1785`)

`rotation_from_direction_roll` computes `finalQuat = alignQuat * rollQuat`.
`alignQuat` maps `dir → X` (world→aligned). Composing the roll on the right
means the roll is applied **about the world X axis before aligning**, so the
combined R no longer maps the material direction onto local X:
`R·dir = A·Rx(θ)·dir ≠ X` unless `dir ∥ X` or `θ = 0`.
Consequence: with a nonzero roll angle and `dir` not along world X, the
stretch is applied along a *wrong* direction. Fix:

    Eigen::Quaternionf finalQuat = rollQuat * alignQuat;   // roll after align

(The bug is invisible when `stretch.y == stretch.z`, because then the roll is
a symmetry of the metric — that's probably why it went unnoticed.)

### 1.3 BUG — background metric ignores the global roll angle

In `compute_scalar_field` (GeneratorLewiner.cpp:328-331) the background
source gets `direction = anisotropyVec` and the stretches, but **`angle` is
never set to `anisotropyAngle`** (stays 0). The uniform path *does* use
`anisotropyAngle` (line 341). So switching from 0 sources to ≥1 source
silently drops the global roll.

### 1.4 BUG (latent) — constructor discards all its arguments

`AnisotropySource(center, direction, angle, stretch, sigma){}`
(Anisotropy.h:24-30) has an **empty body and no member-init list** — every
argument is ignored. Currently all call sites use `make_shared<AnisotropySource>()`
so it's latent, but the first parameterized construction will silently produce
a default source. Also `M` is an `Eigen::Matrix3f` that is **uninitialized
memory** until `update_metric()` is called; `compute_scalar_field` covers this
by calling `update_metric()` on every source, but any other path reading `M`
first gets garbage. Initialize `M = Matrix3f::Identity()` and fix the ctor.

### 1.5 Robustness — `Mblend.inverse()` can produce NaNs

The `safe_sq` clamp was commented out (Anisotropy.h:126-130). The GUI
(`InputFloat`) allows `stretch = 0` (or negative), which makes C singular and
`inverse()` NaN — the whole scalar field becomes NaN. Restore the clamp
(`max(1e-6, s²)`) or clamp stretch in the GUI.

### 1.6 Blending-scheme remarks (for the paper)

* **Partition of unity**: for `ΣW < 1` the background fills the deficit
  (`C = ΣwᵢCᵢ + (1−W)C₀`), for `ΣW ≥ 1` it's a normalized average. This is
  continuous at W = 1 (both branches give `ΣwᵢCᵢ`) but only C⁰ — there is a
  derivative kink. Cleaner and equivalent in spirit:
  `C = (ΣwᵢCᵢ + w_b C₀) / (Σwᵢ + w_b)` with a fixed background weight
  `w_b = 1` (or the existing `backgroundWeight`). Note the header member
  `GeneratorLewiner::backgroundWeight` and its doc comment refer to the *old,
  commented-out* blend and are now dead — remove or rewire.
* **Truncation discontinuity**: the hard cutoff at r > 3σ discards a residual
  weight of exp(−4.5) ≈ 0.011 — a small C⁰ jump in the field at the source
  boundary that the smoothing pass then has to hide. Either subtract the tail
  (`w = exp(−r²/2σ²) − exp(−4.5)`, clamped at 0) or use a compactly supported
  C² kernel (e.g. Wendland).
* **Linear SPD blending swells**: linearly interpolating SPD tensors inflates
  the determinant between differently-oriented sources ("swelling effect").
  Blending C and inverting (≈ harmonic mean of metrics) is a defensible
  choice, but a reviewer may ask about it — the standard alternative is
  log-Euclidean blending: `M(x) = exp(Σ wᵢ log Mᵢ / Σ wᵢ)`. At minimum, state
  the choice explicitly in the paper.

### 1.7 Approximation — Euclidean kNN + metric re-ranking is not exact

In the varied path (GeneratorLewiner.cpp:447-474) the k candidates come from
a *Euclidean* kNN query and are then re-ranked with the local metric. The true
metric-nearest seed can lie outside the Euclidean candidate set; the bound
`k = 3κ²` (κ = global stretch ratio, `choose_candidate_number`) is a heuristic,
not a guarantee, and it also makes cost grow quadratically with κ (κ = 5 →
75 candidates per voxel). This is a reasonable engineering trade-off — say so
in the paper ("candidate set of size 3κ² was found sufficient empirically").

### 1.8 Crash risk — fewer than 3 seeds

`std::partial_sort(cand.begin(), cand.begin()+3, …)` and `cand[2]` /
`neighbors[2]` are UB if fewer than 3 seeds exist. Add an early-out
(`if (seeds.size() < 3) return;` with a logged error).

---

## 2. Porosity vs. container volume

Current formula (GeneratorLewiner.cpp:3473-3474):
`porosity = 1 − V_mesh / V_container` where `V_mesh` is the divergence-theorem
volume of the marching-cubes mesh and `V_container` is the *analytic* (box,
cylinder) or exact mesh volume.

The numerator and denominator live on **different measurement bases**, so the
porosity carries a systematic, resolution-dependent bias — the mesh is always
strictly smaller than the container:

1. `seal_grid_boundaries()` (line 2509) forces the outer 1-voxel shell to air.
   For a box container that coincides with the grid bounds, the scaffold is
   shaved by ≥1 voxel on every face.
2. Marching cubes places the surface by trilinear interpolation inside the
   boundary cells — sub-voxel inset everywhere.
3. Smoothing (any variant) slightly thins struts before extraction.
4. `surfaceMargin = 2.0f` (line 384) is in absolute world units, not voxels —
   at coarse resolutions or scaled scenes the band can interact with boundary
   cells differently.

**Recommended fix**: compute both terms on the *same grid*. You already cache
`containerDistField`; then

    N_domain = #{ x : containerDist(x) ≤ 0 }
    N_solid  = #{ x : scalarField(x) < isoLevel }
    porosity = 1 − N_solid / N_domain

Discretization bias then largely cancels between numerator and denominator and
the estimate converges as resolution increases. Keep the mesh-based value as a
secondary "apparent porosity" if you want it, but report the voxel one.

### 2.1 BUG — `_update_bounding_box()` shrinks the generation domain

`_update_bounding_box()` (line 908) overwrites `bounds[0..5]` — the *grid
definition* — with the **mesh AABB** at the end of every `marching_cubes()`.
Since the mesh AABB is strictly inside the container (points 1-3 above), the
next `compute_scalar_field()` call recomputes `stepX/Y/Z` from the shrunken
bounds (`update_steps`, line 2561) — **the domain erodes on every
regeneration** (change threshold → regenerate → smaller scaffold → …).
Store the mesh AABB only in `aabb` and stop writing into `bounds`.

---

## 3. Scalar-field smoothing options

| | Box 3×3×3 (commented) | Binomial 1-2-1 "Gaussian" | Taubin λ/μ (active) |
|---|---|---|---|
| Frequency response | sinc — negative lobes, can ring | positive, monotone low-pass | band-pass (shrink+inflate) |
| Shrinkage of struts | strongest | moderate | ≈ none (by design) |
| Grid-axis anisotropy | worst (cubic kernel) | mild | mild |
| Effect on porosity/thickness metrics | biases both down | biases both down | approximately neutral |

Taubin is the right default for a scaffold whose *thickness and porosity are
the reported quantities* — the λ/μ pair (0.5 / −0.53) satisfies the stability
condition 0 < λ < −μ, pass-band k_PB = 1/λ + 1/μ ≈ 0.113 (typical range
0.01–0.1, so slightly aggressive but fine). Report (n_iter, λ, μ) in the paper
for reproducibility.

Implementation issues found:

1. **The skip-optimization never fires.** Outside voxels are set to exactly
   `isoLevel + 1.0f` (line 427), and both smoothers test
   `scalarField[idx] > isoLevel + 1.0f` — *strictly* greater. So the entire
   exterior is convolved every iteration: wasted time (the majority of voxels
   for a mesh container) and the exterior plateau slowly erodes. Change to
   `>=` (or set the exterior sentinel to `isoLevel + 1.5f`). Applies to
   `smooth_scalar_field` (line 628) and both Taubin passes (5314, 5343).
2. **Anisotropic voxels**: the 6-neighbor Laplacian `/6` treats stepX/Y/Z as
   equal. If the grid spacing is non-cubic, smoothing is anisotropic in
   physical space. Either enforce cubic voxels or weight neighbors by 1/step².
3. Field-space Taubin is not the mesh-based Taubin of the original paper; the
   non-shrinkage argument transfers only approximately. If a reviewer pushes,
   the clean statement is: "a λ|μ band-pass filter applied to the discrete
   scalar field prior to isosurface extraction". Pipeline order (smooth →
   boundary clamp → island removal → seal) is correct — the clamp restores a
   sharp container boundary after smoothing.
4. Smoothing can change topology (merge/split thin struts) — since you report
   connectivity density, always compute metrics *after* smoothing (currently
   true) and mention it.

---

## 4. OpenMP `collapse(3)` warning

The warning is MSVC **C4849: "OpenMP 'collapse' clause ignored"**. CMake's
`OpenMP::OpenMP_CXX` on MSVC adds `/openmp` = classic OpenMP **2.0** runtime;
`collapse` is an OpenMP 3.0 feature, so the clause is dropped.

Consequences:
* **Correctness: none.** The loop still parallelizes — over the *outer* index
  only.
* **Performance:** irrelevant where the outer extent ≥ thread count (main
  voxel loops, dims ~100+). It *does* hurt in
  `classify_container_narrow_band` (lines 204, 235): the coarse grid outer
  extent is `⌈(N−1)/8⌉+1` ≈ 13–26 for typical dims — fewer chunks than threads
  on a modern CPU, so cores idle in the one phase that calls the expensive
  mesh SDF.

Options:
1. `/openmp:llvm` (VS2019 16.9+; with CMake ≥ 3.30 set
   `set(OpenMP_RUNTIME_MSVC llvm)` before `find_package(OpenMP)`, otherwise
   append the flag manually). Recent VS 2022 versions accept `collapse` under
   the LLVM runtime — verify the warning disappears on your toolset.
2. Toolchain-proof: manually flatten the hot loops —
   `#pragma omp parallel for` over a single `int64 v = 0..total-1` and decode
   `i,j,k` from `v`. Also improves load balance for the narrow-band pass.
3. Do nothing but silence: keep `collapse` for GCC/Clang and accept
   outer-loop-only parallelism on MSVC. If so, at least flatten the two
   coarse-grid loops in `classify_container_narrow_band`.

Note: OpenMP 2.0 also requires *signed* loop indices — all loops already use
`int`, good.

---

## 5. LaTeX pseudocode (algorithm + algpseudocode)

Preamble: `\usepackage{algorithm}` `\usepackage{algpseudocode}` `\usepackage{amsmath,amssymb}`

```latex
\begin{algorithm}[t]
\caption{Anisotropic Voronoi-based scaffold generation}
\label{alg:scaffold}
\begin{algorithmic}[1]
\Require Container $\Omega$ with signed distance $\phi_\Omega$; seeds
  $P=\{p_1,\dots,p_n\}\subset\Omega$; voxel grid $G$ ($N_x\times N_y\times N_z$)
  covering $\Omega$; openness $t\in[0,1]$; thickness map $\tau(\cdot)$;
  background frame $(R_0,S_0)$; anisotropy sources
  $\{(c_i,R_i,S_i,\sigma_i)\}_{i=1}^{m}$; iso-level $\ell$;
  Taubin parameters $(n_s,\lambda,\mu)$, $0<\lambda<-\mu$
\Ensure Watertight triangle mesh $\mathcal{M}$ of the scaffold
\State $C_i \gets R_i^{\top} S_i^{2}\, R_i$ \textbf{for} $i=0,\dots,m$
  \Comment{structure (inverse-metric) tensors}
\State Classify voxels of $G$ as \textsc{Inside} / \textsc{Outside} /
  \textsc{NearSurface} using a coarse Lipschitz bound on $\phi_\Omega$
  \Comment{narrow band; exact SDF only near $\partial\Omega$}
\ForAll{voxels $x \in G$ \textbf{in parallel}}
  \State $\delta \gets \phi_\Omega(x)$
  \If{$\delta > \delta_{\max}$}
     $f(x) \gets \ell + 1$ \Comment{far outside: air} \State \textbf{continue}
  \EndIf
  \State $\tau_x \gets \tau(|\delta|)$ \Comment{local wall thickness (graded or uniform)}
  \State $w_i \gets \exp\!\big(-\lVert x-c_i\rVert^{2}/2\sigma_i^{2}\big)$,\quad
         $W \gets \textstyle\sum_{i=1}^{m} w_i$
  \If{$W \ge 1$}
     \State $C(x) \gets \tfrac{1}{W}\sum_i w_i\, C_i$
  \Else
     \State $C(x) \gets \sum_i w_i\, C_i + (1-W)\, C_0$
     \Comment{background fills the deficit}
  \EndIf
  \State $M(x) \gets C(x)^{-1}$ \Comment{local Riemannian metric}
  \State $\mathcal{K} \gets$ $k$ Euclidean nearest seeds of $x$
    (kd-tree), $k = \lceil 3\kappa^2 \rceil$,
    $\kappa = s_{\max}/s_{\min}$
  \State Re-rank $\mathcal{K}$ by
    $d_M(x,p)^2 = (x-p)^{\top} M(x)\,(x-p)$;\;
    keep the three smallest $d_1 \le d_2 \le d_3$
  \State $\nabla d_j \gets M(x)\,(x-p_j)\,/\,d_j$, \quad $j=1,2,3$
  \State $v \gets (1-t)\,(d_2-d_1) \;+\; t\,(d_3-d_1)$
    \Comment{$t{=}0$: walls (foam), $t{=}1$: struts (lattice)}
  \State $g \gets (1-t)\,(\nabla d_2-\nabla d_1) + t\,(\nabla d_3-\nabla d_1)$
  \State $f(x) \gets 2\,v\,/\max(\lVert g\rVert,\varepsilon)
    \;-\; \tau_x \;+\; \ell$
    \Comment{gradient-normalised, thickness-shifted}
\EndFor
\For{$s = 1,\dots,n_s$} \Comment{Taubin $\lambda|\mu$ band-pass smoothing}
  \State $f \gets f + \lambda\,\Delta_6 f$;\qquad
         $f \gets f + \mu\,\Delta_6 f$
         \Comment{$\Delta_6$: 6-neighbour discrete Laplacian}
\EndFor
\State $f(x) \gets \max\!\big(f(x),\; \phi_\Omega(x) + \ell\big)$
  \Comment{Boolean intersection with $\Omega$}
\State Remove disconnected solid components of $\{f < \ell\}$ (BFS, keep largest);
  set the outer grid shell to air \Comment{watertightness}
\State $\mathcal{M} \gets \textsc{MarchingCubes33}(f,\ell)$
  \Comment{topologically consistent MC~\cite{Lewiner2003}}
\State \Return $\mathcal{M}$
\end{algorithmic}
\end{algorithm}
```

Notes for the manuscript:
* Cite Lewiner et al. 2003 (MC33) for the extraction step and Taubin 1995 for
  the λ|μ filter.
* State the uniform-anisotropy special case as an optimization: when m = 0 the
  metric is constant, so instead of per-voxel re-ranking the seeds are warped
  once by $S_0^{-1}R_0$ and a standard Euclidean kd-tree query is used —
  mathematically identical, much cheaper.
* The factor 2 in the normalization makes $f$ approximate the full local
  wall/strut thickness rather than the half-width; if you drop it in the code,
  drop it here too.

---

## 6. Issue summary (ranked)

**Status update (2026-07-09, second pass):**
- #2 fixed — `_update_bounding_box` no longer writes `bounds`; `get_bounds()`
  now returns the mesh AABB (what all external callers wanted); the STL-load
  constructor sets `bounds` explicitly. This also fixes latent grid/AABB
  misalignment in `extract_from_ROI` and `estimate_anisotropy` ROI indexing.
- #3 fixed — primary `porosity` is now voxel-based (computed in
  `compute_scalar_field` on the shared grid); mesh-based value kept as
  `porosityMesh` (legacy), shown in the GUI and exported as "Mesh Porosity".
- #4 fixed by user (background.angle now set).
- #5 fixed by user (`safe_sq` restored in both metric builders).
- #6 fixed — skip tests are `>=`; Taubin Laplacian now weights axes by
  1/step² (isotropic in world space for non-cubic voxels; identical to the
  old /6 average for cubic voxels).
- #7 fixed — ctor now uses its arguments and calls `update_metric()` only.
  `update_model()` was removed from the ctor again: it creates OpenGL buffers
  and would crash the CLI profiler / the .scaf loader thread (both construct
  sources without a GL context on the current thread).
- #8 fixed by user (seed-count guard; logger pointer now null-checked).
- #1 (roll composition) — user reports the direction looks correct in
  practice; decisive test proposed (dir = +Y, roll = 90°, stretch = (3,1,1):
  current code predicts elongation along Z, correct code along Y).
- #9 (OpenMP collapse) — deferred by user decision.
- #10/#11/#12 — open (kernel truncation, dead `backgroundWeight`, duplicated
  metric builders).


| # | Severity | Where | Issue |
|---|----------|-------|-------|
| 1 | **High** | `Utils.cpp:1785` | Roll quaternion composed on the wrong side — roll about world X instead of the material direction; wrong metric whenever roll ≠ 0 and dir ∦ X and stretch.y ≠ stretch.z |
| 2 | **High** | `GeneratorLewiner.cpp:908` (`_update_bounding_box`) | Grid `bounds` overwritten with mesh AABB → generation domain shrinks on every regeneration |
| 3 | **High** | `estimate_metrics` :3473 | Porosity mixes mesh volume with analytic container volume — systematic overestimate; compute both on the voxel grid |
| 4 | Medium | `compute_scalar_field` :328 | Background anisotropy ignores `anisotropyAngle` when sources exist |
| 5 | Medium | `Anisotropy.h:254` | `inverse()` of possibly singular C (stretch = 0 allowed in GUI); `safe_sq` clamp removed |
| 6 | Medium | Taubin :5314/:5343, Gaussian :628 | Skip test `> isoLevel + 1.0f` never true for exterior voxels (they equal it exactly) — whole exterior smoothed every iteration |
| 7 | Medium | `Anisotropy.h:24` | Constructor silently discards all arguments; `M` uninitialized until `update_metric()` |
| 8 | Medium | varied path :461 | `partial_sort(...begin()+3...)` UB when < 3 seeds/candidates |
| 9 | Low | CMake / MSVC | `collapse(3)` ignored (C4849) — perf-only; hurts mainly the coarse narrow-band loops (outer extent ~13–26) |
| 10 | Low | `Anisotropy.h:234` | Hard 3σ cutoff leaves a 0.011 weight discontinuity in the blended field |
| 11 | Low | header :352 | `backgroundWeight` member + doc comment are dead (refer to the removed blend variant) |
| 12 | Low | `Anisotropy.h` | `create_metric()` free function duplicates `update_metric()` — one will drift |
| 13 | Low | :384 | `surfaceMargin = 2.0f` in absolute units, not scaled to voxel size |
| 14 | Info | :5330 | 6-neighbor Laplacian assumes cubic voxels (stepX=stepY=stepZ) |
