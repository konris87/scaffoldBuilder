# Third-party licenses

ScaffoldBuilder bundles the following third-party libraries under `lib/`. Each is
distributed under its own license, reproduced or referenced below. All are
compatible with ScaffoldBuilder's GPL-3.0 licensing.

> **Maintainer note (pre-release):** two entries below are marked ⚠️ — their
> license *text* is not currently vendored in the tree (only headers / generated
> code are). Before a binary or source release, add the missing license files so
> the distribution is self-contained, as GPLv3 §5/§6 and each permissive license
> require.

| Library | Version | License | License text in tree |
|---|---|---|---|
| Dear ImGui | (vendored) | MIT | `lib/imgui/LICENSE.txt` |
| GLFW | (vendored) | Zlib/libpng | `lib/glfw/LICENSE.md` |
| GLM | 1.0.1 | MIT / The Happy Bunny License | `lib/glm-1.0.1/copying.txt` |
| ImGuiFileDialog | (vendored) | MIT | `lib/ImGuiFileDialog/LICENSE` |
| Eigen | (vendored, headers) | MPL-2.0 (core) | ⚠️ not vendored — see below |
| glad | generated | Public domain / MIT (KHR headers) | ⚠️ not vendored — see below |

---

## Dear ImGui — MIT
Omar Cornut and contributors. Full text: [`lib/imgui/LICENSE.txt`](lib/imgui/LICENSE.txt).

## GLFW — Zlib/libpng
Copyright © Marcus Geelnard and Camilla Löwy.
Full text: [`lib/glfw/LICENSE.md`](lib/glfw/LICENSE.md).

## GLM — MIT / The Happy Bunny License
Copyright © G-Truc Creation. ScaffoldBuilder uses it under the **MIT** option.
Full text: [`lib/glm-1.0.1/copying.txt`](lib/glm-1.0.1/copying.txt).

## ImGuiFileDialog — MIT
Copyright © Stephane Cuillerdier (aiekick).
Full text: [`lib/ImGuiFileDialog/LICENSE`](lib/ImGuiFileDialog/LICENSE). (Its bundled
`dirent` port carries its own MIT license at `lib/ImGuiFileDialog/dirent/LICENSE`.)

## Eigen — MPL-2.0  ⚠️
Only the `Eigen/` headers are vendored; the upstream top-level license files are
not. The parts ScaffoldBuilder uses (`Eigen/Dense`, geometry, decompositions) are
**MPL-2.0**. Some optional modules Eigen ships are LGPL/BSD, but none of those are
compiled here.
- **Action:** add Eigen's `COPYING.MPL2` (from the Eigen 3.4 source release) as
  `lib/Eigen/COPYING.MPL2`, or record the pinned Eigen version and upstream URL
  (https://eigen.tuxfamily.org / https://gitlab.com/libeigen/eigen).

## glad — generated loader  ⚠️
`lib/glad` is output from the glad generator (https://glad.dav1d.de,
`profile=core, api=gl 4.6, loader=on`). The generated loader code is placed in the
**public domain** by the generator; the bundled Khronos `KHR/khrplatform.h` and GL
headers are under the **MIT / Apache-2.0** Khronos license.
- **Action:** add the Khronos header license (SPDX `MIT` header block already
  present in the KHR headers) and a short note that the loader is generated /
  public domain.

---

*Add ScaffoldBuilder's own `LICENSE` (GPL-3.0) at the repository root — it is
referenced by `README.md` but not yet present.*
