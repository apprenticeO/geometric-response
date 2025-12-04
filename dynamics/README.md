# Scripts_Dyanmics_triad

Minimal computational companions for “From Structural Triad to Response Theory”.

- 01_xxz_pitriad_ed.py
  - Closed XXZ chain (small N) via exact diagonalization
  - Computes time series of Π_A(t) for a block A, its autocorrelation C_Π(t), Green–Kubo time τ_G, and χ_Π(ω) by FFT
  - Outputs CSVs under the same folder

Conventions:
- ħ = 1, pure global evolution ⇒ I(A:Ā)=2 S_A
- √F_Q(ρ_A;H_A)=2 Δ_{ρ_A} H_A (Mandelstam–Tamm/SLD)
- Default parameters are defined at top of each script; run with `python3` (no flags)


