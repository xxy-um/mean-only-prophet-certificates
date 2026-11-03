# Verification records

These records document checks of the exact certificate data. The exact integer
numerators and denominators determine acceptance; JSON decimal numbers are
rounded displays.

- `RELEASE_FULL_REPLAY.log` records a complete replay of all 5,120 intervals.
- `GENERATION_FROM_RELEASE.json` records reconstruction of all 5,120 proposal
  files from the frozen seeds, with byte-for-byte comparison.
- `critical_ubsan.jsonl` and `CRITICAL_REPLAY_REPORT.json` record four selected
  arithmetic checks with the undefined-behavior sanitizer enabled. Their source
  paths are relative to `ps/`.
- The other records retain the original independent input, scope, target, and
  coverage checks. File counts and timings describe their recorded runs; the
  current release file list is in `../certificate.json`.

The supplementary proof in `../Verification.pdf` supplies the implementation
and interval-arithmetic justifications. Numerical replay and checksum checks
complement that proof.
