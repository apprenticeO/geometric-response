This repository formalizes results from:

> **From Microscopic Structural Interdependence to Causal Dispersive Geometry**  
> Ovidiu Tataru, Zenodo, 2026  
> https://doi.org/10.5281/zenodo.18832967

This paper builds on the framework introduced in two companion papers 
published in the *International Journal of Quantum Foundations*, vol. 11 (2025):

> **A Quantum Structural Triad: Fluctuations, Entropy, and Correlations 
> as Interdependent Primitives**  
> Ovidiu Tataru, IJQF 11, 800–852 (2025)  
> https://ijqf.org/wp-content/uploads/2025/11/IJQF2025v11n4p20r.pdf  
> *Includes companion Rocq/Coq formalization — 
> github.com/apprenticeO/Structural-Triad*

> **Geometric Response and the Frequency-Dependent Gravitational Coupling**  
> Ovidiu Tataru, IJQF 11 (2025)  
> https://ijqf.org/wp-content/uploads/2025/11/IJQF2025v11n4p21.pdf

The paper develops a Green–Kubo → Kramers–Kronig chain from the microscopic
structural observable \(\Pi(t)=\sqrt{F_Q}\,S\,I\) to a causal, frequency–dependent
geometric response \(G_{\mathrm{eff}}(\omega)=G_0\,\chi(\omega)\).  
This repository packages:

- a complete Lean 4 formalization of the analytic and conservation theorems, and  
- the minimal Python tooling and sample data needed to reproduce the GK/KK
  diagnostics and Debye fits used in the manuscript.

---

## Contents

- `Geometric_Response_Lean/`  
  Lean 4 project `GeometricResponseLean`:
  - Green–Kubo to Kramers–Kronig necessity (GK→KK small‑\(\omega\) slope),
  - Debye single–pole sufficiency and passivity,
  - flat–space conservation for polynomial \(G(\Box)\),
  - convolution‑kernel conservation (CTP viewpoint),
  - Wiener–Khinchin \(S(0)=2\tau_G\),
  - Paley–Zygmund positivity and \(\tau_G>0\),
  - plus two explicit physics axioms (short‑time variance bound and \(L^1\)
    autocorrelation under H1–H3).

- `scripts/`  
  Geometry–response pipeline and diagnostics:
  - `generate_ou.py`, `generate_lindblad_qubit.py`: construct OU and Lindblad
    benchmark time series \(\Pi(t)\),
  - `compute_gk_kk.py`, `fit_debye.py`, `write_debye_chi.py`, `write_scales.py`:
    Green–Kubo integrals, KK slopes, Debye fits, and derived scales
    \(\{\tau_G,\omega_G,\lambda_G,m_G\}\),
  - `chain_two_qubit_pipeline.py`, `run_two_qubit_demo.py`,
    `test_two_qubit_lwss_passivity.py`, `verify_transfer_identity.py`:
    two–qubit structural tests and passivity checks,
  - `build_tensor_transfer_dataset.py`, `synthesize_response_from_rinfo.py`,
    `export_rmunu_diagnostic.py`, `calibrate_G0.py`:
    tensor response datasets and \(G_0\) calibration,
  - `utils.py`, `QuantumTriadVisualizer.tsx`: shared helpers and visualization.

- `dynamics/`  
  Minimal triad‑dynamics scripts that generate representative \(\Pi(t)\)
  trajectories used as structural input for the geometry pipeline:
  TEBD/XXZ examples, active–epoch density tests, and Paley–Zygmund validation.

- `tests/`  
  Stand‑alone diagnostics that mirror the “assumptions & tests” table in the
  paper (LWSS / stationarity, flatness and \(\tau_{\mathrm{PSD}}\),
  finite‑record bias, multi‑window \(\tau_G\), and violation cases).

- `data/`  
  Small, curated CSV/PNG snapshots for:
  - `baseline/`: OU baseline GK/KK checks and Debye fit,
  - `lindblad/`: Lindblad‑qubit GK/KK checks,
  - `violation/`: violation summary and plots.
  Large raw time–series files are **not** included; they can be regenerated
  by running the scripts.

---

## How This Relates to the Paper

- Theorems about the GK→KK link, Debye sufficiency, conservation, and
  Wiener–Khinchin are **fully machine‑verified** in Lean 4
  (see `Geometric_Response_Lean/PROVEN_VS_PAPER.md` for a line‑by‑line map
  between LaTeX statements and Lean theorems).
- The Python scripts implement exactly the GK/KK estimator pipeline,
  Debye fits, and diagnostic plots described in Appendices A–B of the paper.
- The `data/` directory contains the minimal numerical artifacts needed to
  reproduce the figures and tables referenced in the appendices without
  shipping multi‑hundred‑MB raw traces.

---

## Third‑Party Software and Dependencies

This repository **depends on** but does not redistribute:

- **Lean 4** and **Mathlib** – used by `GeometricResponseLean`.  
  They are available under their respective open‑source licenses from  
  `https://github.com/leanprover/lean4` and `https://github.com/leanprover-community/mathlib4`.

- **Python** and common scientific packages (e.g. `numpy`, `scipy`, `matplotlib`,
  `pandas`) – used by the scripts in `scripts/`, `dynamics/`, and `tests/`.

All such third‑party projects remain under their own licenses. This repository
only includes *glue code and models* that call into those libraries; it does
not modify or redistribute them.

---

## Copyright and Disclaimer

Copyright © 2025 Ovidiu Tataru.

Permission is granted to read and use the code and proofs in this repository
for research and educational purposes. Redistribution or derivative works
should preserve attribution to the original author and clearly indicate any
modifications.

This repository is provided **“as is”**, without any warranty of correctness
or fitness for a particular purpose. The Lean formalization is meant to check
the internal consistency of the mathematical framework; it does **not**
constitute an endorsement or validation of the physical modelling assumptions.
All third‑party tools and libraries are used under their own licenses and
remain the responsibility of their respective authors and maintainers.


