# Independent row checks

These files audit selected rows of the **original Appendix A.1 LP**. They are additional checks; the global feasibility argument is the backward construction and safe coverage of every original rate leaf described in [`../../docs/verification.pdf`](../../docs/verification.pdf).

`reference_rows.cpp` independently scans all 131,093 original rate rows for each requested state using long-double arithmetic and the frozen coefficient approximations. It is a diagnostic scan, not a rigorous interval proof. `check_exact_row.py` instead reconstructs the exact coefficients from rational inputs and encloses exponentials, logarithms, and the final row using rational intervals. `check_selected_exact.py` applies that interval check to the three closest rows in a scan and to the distinguished state `(1200,0)`; it handles the zero-rate successor identity exactly.

Run the following commands from the package root after the replay has produced the referenced checkpoint:

```sh
g++ -O2 -std=c++17 -I data audit/reference_rows.cpp -o audit/reference_rows
./audit/reference_rows data/grids/grid1200.txt runs/run/layer_320.next 320 \
    < audit/potential16_final320_states.txt > runs/run/rescan_320.jsonl
python3 audit/check_selected_exact.py runs/run/rescan_320.jsonl \
    runs/run/layer_320.next runs/run/exact_rescan_320.json
```

Run the Python row checkers without `-O` or `PYTHONOPTIMIZE`; they reject disabled assertions.

The scanner input is one `i j post_integer` triple per line. The `.next` file is the already fixed **next-layer** table, not the current `.table` file. The stored state triples belong to the recorded run; use matching checkpoint data.

For a single row, supply a JSON object with `k`, `i`, `j`, the zero-based rate `leaf`, `post_integer`, and `table_path`:

```sh
python3 audit/check_exact_row.py audit/candidate_320_example.json
```

A relative `table_path` is resolved against the candidate JSON's directory. The included example points to `runs/run/layer_320.next`. The `*_results.jsonl` and `*_exact.json` files record completed checks at selected time layers. Historical file and command paths have been made relative to their original workspace root; they are descriptive labels, not paths required by the portable scripts. Only the workspace prefix was removed: the recorded arithmetic values, timestamps, and result hashes are unchanged. Pass the locally reconstructed checkpoint path to `check_selected_exact.py` as shown above.
