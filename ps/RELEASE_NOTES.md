# Release notes

## Certified result

| Quantity | Value |
| --- | --- |
| Certified minimum | `49839735470020217/72057594037927936` |
| Decimal value | `0.69166527325054429076...` |
| Acceptance target | `1729/2500 = 0.6916` |
| Exact margin over the target | `2939645868285289/45035996273704960000` |
| Threshold height | `567246577807/800000000000 = 0.70905822225875` |
| Threshold intervals `N` | `160` |
| Stopping-time intervals `N_t` | `5120` |
| Probability-averaging intervals `N_p` | `10240` |

The complete supplementary argument is consolidated in `Verification.tex` and `Verification.pdf`. The source is a single self-contained TeX file.

The numerical parameters, proposal files, arithmetic records, C++ verifier, and deterministic assembly rule are unchanged. The Python launchers validate inputs explicitly, including under optimized Python execution, and the replay command supports repository paths containing spaces. These maintenance changes do not change the certified mathematical instance or acceptance target.

## Manuscript consistency

Section 4.4 and Appendix B.2 of the accompanying manuscript use `N=160` positive nondecreasing drop rates. Appendix B.2 uses `N_t=5120` and `N_p=10240`, consistently with this certificate.

The certified competitive ratio remains `0.6916`. The longer decimal denotes the exact certificate minimum displayed above; it is not a claim that the threshold family or the number of threshold intervals is optimal.

## Reproducibility scope

`replay.py --jobs 8` recompiles the verifier and reproduces all 5,120 directed-arithmetic records. `generation/rebuild_proposals.py` reconstructs all final proposal files from the supplied frozen candidate lists and rational selection schedule. It does not restart the exploratory optimization search from an unspecified initialization. Verification of the feasible certificate does not require that exploratory search.

The archive is complete for local verification.
