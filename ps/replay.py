#!/usr/bin/env python3
"""Audit and replay the fixed, article-only 160-segment FRLP certificate."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
TARGET = Fraction(1729, 2500)
KEYS = ('N', 'lo', 'width', 'lower_numerator', 'denominator',
        'constant_lo', 'constant_hi', 'min_suffix', 'min_suffix_index',
        'max_gradient_width')


def require(condition, message):
    """Validate certificate data even when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def audit():
    manifest = json.loads((ROOT / 'certificate.json').read_text())
    require(manifest['target_numerator'] == TARGET.numerator, 'target numerator mismatch')
    require(manifest['target_denominator'] == TARGET.denominator, 'target denominator mismatch')
    for rel, digest in manifest['sha256'].items():
        require(hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest,
                f'hash mismatch: {rel}')
    parameters = json.loads((ROOT / 'parameters.json').read_text())
    intervals = json.loads((ROOT / 'intervals.json').read_text())
    records = json.loads((ROOT / 'records.json').read_text())
    require(len(intervals) == len(records) == 5120, 'expected 5120 intervals and records')
    expected = None
    for source in sorted({x['source'] for x in intervals}):
        require(not Path(source).is_absolute() and '..' not in Path(source).parts
                and not any(c.isspace() for c in source),
                f'expected a relative proposal path without whitespace: {source}')
        rows = (ROOT / source).read_text().splitlines()
        require(rows[0] == 'MEAN_ONLY_AFFINE_V1', f'proposal format: {source}')
        coarse, k = map(int, rows[1].split())
        require(k == 160 and coarse > 0 and coarse % 160 == 0 and 10240 % coarse == 0,
                f'proposal grid dimensions: {source}')
        params = tuple(map(int, rows[2].split()))
        require(len(params) == 640, f'common parameter dimensions: {source}')
        if expected is None:
            expected = params
        require(params == expected, f'different common parameters: {source}')
        b, m, kap, eta = (params[j:j+160] for j in range(0, 640, 160))
        require(all(v > 0 for v in b) and all(x < y for x, y in zip(b, b[1:])),
                f'positive increasing drop rates: {source}')
        require(sum(b) <= 160 * 10**10, f'threshold height exceeds one: {source}')
        require(all(0 <= x <= y for x, y in zip(m, b)), f'allocation parameters: {source}')
        require(all(x >= 0 for x in kap) and all(x == 0 for x in eta),
                f'allocation or auxiliary parameters: {source}')
        tail = list(map(int, rows[3].split()))
        require(len(tail) == coarse and tail == sorted(tail), f'auxiliary tail dimensions: {source}')
        require(all(0 <= x <= 2**20 for x in tail), f'auxiliary tail range: {source}')
        count = int(rows[4])
        require(len(rows) == count + 5, f'plane count mismatch: {source}')
        mass = [0] * (3 * coarse)
        for line in rows[5:]:
            v = line.split()
            owner, weight = int(v[0]), int(v[1])
            require(0 <= owner < 3 * coarse and 0 <= weight <= 2**32,
                    f'plane owner or weight: {source}')
            mass[owner] += weight
            if v[2] == 'log':
                require(len(v) == 6 and owner < 2 * coarse and int(v[3]) == 1,
                        f'logarithmic plane format: {source}')
            elif v[2] == 'C':
                require(len(v) == 4 and owner >= 2 * coarse and 0 <= int(v[3]) <= 10**10,
                        f'survival plane format: {source}')
            else:
                raise ValueError(f'unsupported plane kind: {v[2]}')
        require(all(x <= 2**32 for x in mass), f'plane mass budget exceeded: {source}')
    require(list(expected[:160]) == parameters['b_numerators'], 'drop-rate parameters mismatch')
    require(list(expected[160:320]) == parameters['m_numerators'], 'allocation parameters mismatch')
    require(list(expected[320:480]) == parameters['kappa_numerators'], 'kappa parameters mismatch')
    require(Fraction(sum(expected[:160]), 160 * 10**10) == Fraction(parameters['theta_0']),
            'threshold height mismatch')
    for j, (x, y) in enumerate(zip(intervals, records)):
        require(x['N'] == 10240 and x['lo'] == 2*j and x['width'] == 2,
                f'incomplete or inconsistent interval coverage: {j}')
        require(x['lo'] // 64 == (x['lo'] + 1) // 64, f'interval crosses a threshold knot: {j}')
        require(all(x[k] == y[k] for k in ('source', 'N', 'lo', 'width')),
                f'input and record assignment mismatch: {j}')
        require(y['status'] == 'PASS', f'saved record did not pass: {j}')
        require(Fraction(int(y['lower_numerator']), int(y['denominator'])) >= TARGET,
                f'saved record below target: {j}')
    minimum = min(Fraction(int(y['lower_numerator']), int(y['denominator'])) for y in records)
    require(minimum == Fraction(int(manifest['minimum_numerator']), int(manifest['minimum_denominator'])),
            'manifest minimum mismatch')
    print(f'AUDIT PASS: target={TARGET} ({float(TARGET)}), all 5120 intervals, common parameters, plane budgets, hashes and coverage.', flush=True)
    return intervals, records, minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--critical-only', action='store_true')
    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error('--jobs must be positive')
    intervals, records, minimum = audit()
    if args.audit_only:
        return
    ids = list(range(len(intervals)))
    if args.critical_only:
        worst = min(ids, key=lambda j: Fraction(int(records[j]['lower_numerator']), int(records[j]['denominator'])))
        ids = sorted({0, worst, len(intervals)//2, len(intervals)-1})
    with tempfile.TemporaryDirectory(prefix='ps160-article-replay-') as work:
        work = Path(work)
        exe = work / 'verify'
        subprocess.run(['g++', '-O3', '-std=c++17', str(ROOT / 'sources/strict160_verify.cpp'), '-o', str(exe)], check=True)
        def replay_shard(pair):
            k, selected = pair
            tasks, output = work / f'tasks{k}', work / f'output{k}.jsonl'
            # The fixed proposal names contain no whitespace. Resolve them from
            # ROOT so a space in a parent directory never enters the task format.
            tasks.write_text('\n'.join(f"{intervals[j]['source']} 10240 {intervals[j]['lo']} 2" for j in selected) + '\n')
            with output.open('w') as stream:
                subprocess.run([str(exe), str(tasks), str(TARGET.numerator), str(TARGET.denominator), 'RECTANGLE_V1'],
                               cwd=ROOT, stdout=stream, check=True)
            got = [json.loads(line) for line in output.read_text().splitlines()]
            require(len(got) == len(selected), f'arithmetic replay record count: shard {k}')
            for j, row in zip(selected, got):
                require(row['status'] == 'PASS', f'arithmetic replay did not pass: {j}')
                for key in KEYS:
                    require(row[key] == records[j][key], f'arithmetic replay mismatch: interval {j}, {key}')
            print(f'REPLAY PASS: shard {k}, {len(got)} intervals.', flush=True)
            return len(got)
        jobs = min(args.jobs, len(ids))
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            count = sum(pool.map(replay_shard, [(j, ids[j::jobs]) for j in range(jobs)]))
    print(json.dumps({'status': 'VERIFIED', 'target_numerator': TARGET.numerator,
                      'target_denominator': TARGET.denominator, 'target_met': minimum >= TARGET,
                      'exact_target_margin': str(minimum-TARGET),
                      'complete_input_and_coverage_audit': True,
                      'all_5120_arithmetic_checks_replayed': not args.critical_only,
                      'arithmetic_checks_replayed': count,
                      'minimum_numerator': str(minimum.numerator),
                      'minimum_denominator': str(minimum.denominator),
                      'minimum_decimal': float(minimum)}, indent=2))


if __name__ == '__main__':
    main()
