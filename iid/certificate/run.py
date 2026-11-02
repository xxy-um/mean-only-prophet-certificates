#!/usr/bin/env python3
"""Compile and replay the literal Appendix A.1 FRLP backward evaluator."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument('--output', type=Path, default=ROOT / 'runs' / 'run')
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--skip-input-validation', action='store_true',
                        help='reuse a previously completed exact input check')
    parser.add_argument('--slack-bits', type=int, default=28)
    args = parser.parse_args()
    if args.threads < 1 or not 20 <= args.slack_bits <= 40:
        parser.error('threads must be positive and slack-bits must lie in [20,40]')
    source = ROOT / 'verify.cpp'
    if not source.exists():
        parser.error('the frozen verify.cpp file is missing')
    output = args.output.resolve()
    signature = {
        'source_sha256': sha256(source),
        'data_sha256': {
            str(p.relative_to(ROOT)): sha256(p)
            for p in sorted(DATA.rglob('*')) if p.is_file()
        },
        'slack_bits': args.slack_bits,
        'potential_bits': 16,
        'N_t': 640, 'N_y': 1200, 'N_q': 131093, 'slopes': 40,
    }
    metadata = output / 'RUN_INPUTS.json'
    if args.resume:
        if not metadata.exists() or not (output / 'latest.k').exists():
            parser.error('resume requires RUN_INPUTS.json and a completed checkpoint')
        if json.loads(metadata.read_text()) != signature:
            parser.error('resume source, inputs, or numerical parameters differ')
    else:
        if output.exists() and any(output.iterdir()):
            parser.error('output directory is not empty; use --resume or another --output')
        output.mkdir(parents=True, exist_ok=True)
        metadata.write_text(json.dumps(signature, indent=2) + '\n')

    if not args.skip_input_validation:
        validation = subprocess.run(
            [sys.executable, str(ROOT / 'validate_inputs.py')],
            text=True, capture_output=True, check=True)
        (output / 'INPUT_VALIDATION.json').write_text(validation.stdout)
        print(validation.stdout, end='')

    binary = output / 'verify'
    compile_command = [
        args.compiler, '-O3', '-std=c++17', '-fopenmp', '-ffp-contract=off',
        '-I', str(DATA), '-DCOEFFICIENT_FILE="coeff_N640_Q131093.h"',
        str(source), '-o', str(binary),
    ]
    with (output / 'compile.log').open('w') as stream:
        subprocess.run(compile_command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    command = [str(binary), str(args.threads), str(DATA / 'grids' / 'grid1200.txt'),
               str(output), str(args.slack_bits), 'resume' if args.resume else 'fresh', '0']
    (output / 'COMMAND.json').write_text(json.dumps({
        'compile': compile_command, 'replay': command,
    }, indent=2) + '\n')
    started = time.monotonic()
    # Logs remain available even if the replay is interrupted; completed
    # checkpoints can be resumed with the same inputs and source.
    with (output / 'replay.stderr').open('a' if args.resume else 'w') as stream:
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=stream)
    (output / 'replay.stdout').write_text(result.stdout)
    if result.returncode:
        raise SystemExit(f'evaluator exited with status {result.returncode}; see replay.stderr')
    rows = [line for line in result.stdout.splitlines() if line.strip()]
    if not rows:
        raise SystemExit('evaluator returned no final JSON result')
    record = json.loads(rows[-1])
    if record.get('last_k') != 0:
        raise SystemExit('the evaluator did not finish the backward run at k=0')
    numerator, denominator = int(record['numerator']), int(record['denominator'])
    if denominator != 2**40:
        raise SystemExit('unexpected result scale')
    record.update({
        'target_met': 10000 * numerator >= 7419 * denominator,
        'exact_target_margin': 10000 * numerator - 7419 * denominator,
        'launcher_seconds': time.monotonic() - started,
    })
    (output / 'RESULT.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))
    if not record['target_met']:
        raise SystemExit('the completed feasible run did not reach 0.7419')


if __name__ == '__main__':
    main()
