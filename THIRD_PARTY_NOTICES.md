# Third-Party Notices

This repository uses the following third‑party tools and libraries.  
This file provides attribution and links to upstream licenses and documentation.
No third‑party source code is modified here; external tools are consumed as
dependencies via their standard distribution channels.

> **Note.** You are not redistributing the Lean/Mathlib toolchain or Python
> packages in this repository; they are installed separately by users.

---

## Lean 4 and Mathlib

- **Lean 4** (theorem prover and programming language)  
  - Upstream: `https://github.com/leanprover/lean4`  
  - License: see `LICENSE` in the Lean repository.
  - Usage: toolchain for building and checking `Geometric_Response_Lean/GeometricResponseLean`.

- **Mathlib4** (Lean mathematical library)  
  - Upstream: `https://github.com/leanprover-community/mathlib4`  
  - License: see `LICENSE` in the Mathlib repository.
  - Usage: measure theory, integration, complex analysis, and spectral tools
    used throughout the Lean development.

Any local cache files or build artifacts created by Lean/Mathlib are generated
at build time and are not part of this repository’s source tree.

---

## Python and Scientific Libraries

The Python scripts under `scripts/`, `dynamics/`, and `tests/` rely on the
following libraries (directly or transitively). Please refer to each project’s
upstream license for current terms.

- **Python** – `https://www.python.org/`
- **NumPy** – `https://numpy.org/doc/stable/license.html`
- **SciPy** – `https://scipy.org/scipylib/license.html`
- **Matplotlib** – `https://matplotlib.org/stable/users/project/license.html`
- **Pandas** – `https://pandas.pydata.org/about/licenses/`

These libraries are typically installed via `pip` or a distribution package
manager. They are not vendored in this repository.

---

## Script‑Specific Dependency Notes (Informative)

The following files use one or more of the libraries listed above. This list is
informational and not exhaustive of transitive dependencies.

- `scripts/generate_ou.py`, `scripts/generate_lindblad_qubit.py`  
  - Uses: NumPy, SciPy (random processes, ODEs), Pandas (CSV I/O), Matplotlib (optional plots).

- `scripts/compute_gk_kk.py`, `scripts/fit_debye.py`,
  `scripts/write_debye_chi.py`, `scripts/write_scales.py`  
  - Uses: NumPy, SciPy, Pandas, Matplotlib (for diagnostic plots).

- `scripts/chain_two_qubit_pipeline.py`, `scripts/run_two_qubit_demo.py`,
  `scripts/test_two_qubit_lwss_passivity.py`, `scripts/verify_transfer_identity.py`  
  - Uses: NumPy, SciPy, Matplotlib, Pandas.

- `scripts/build_tensor_transfer_dataset.py`,
  `scripts/synthesize_response_from_rinfo.py`,
  `scripts/export_rmunu_diagnostic.py`,
  `scripts/calibrate_G0.py`  
  - Uses: NumPy, SciPy, Pandas, Matplotlib.

- `dynamics/*`  
  - Uses: NumPy and Matplotlib (for TEBD/XXZ‑style structural triad diagnostics).

- `tests/*`  
  - Uses: NumPy, SciPy, Pandas, Matplotlib (for LWSS/flatness/bias diagnostics).

---

## Questions

For questions about third‑party notices or to request updates, please contact:

- **Email**: `ovdttr@gmail.com`


