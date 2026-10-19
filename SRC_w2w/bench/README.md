# Performance redesign benchmark ledger

One row per (feature, case, metric). `feature=baseline` rows are the
pre-redesign anchors (development-tree commit fd191ba, not in this
repository's history), measured on Rh3Si4 (logs and design notes are in the
development tree).

Cases:
- `rh3si4_soc`: 352-band window, nk=120, nntot=6, 160 spinor candidates,
  nw=120 (tests/Rh3Si4/w2w_hybrid) — the scaling target.
- `fe_nosoc`:   15-band window, nk=512, 9 proj (tests/w2w_runs/Fe) — the
  fast parity case.

Parity rule: Omega_I to 1e-6 relative and sigma_min to same order of
magnitude vs the baseline anchors (never bitwise once threaded paths are on).
Each increment appends rows tagged with its feature name (A1_telemetry,
A2_colored_omp, A3_fused_monitor, A4_zheevr, A5_mixing, B1_binary_io,
C1_alm_cache, C2_hybr_mode, ...).
