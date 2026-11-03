#!/usr/bin/env python3
"""Verify all 5,120 stopping-time intervals of the fixed PS certificate."""
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
EXPECTED_MINIMUM = Fraction(49839735470020217, 72057594037927936)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit():
    manifest = json.loads((ROOT / 'certificate.json').read_text())
    require(manifest['target_numerator'] == TARGET.numerator, 'target numerator mismatch')
    require(manifest['target_denominator'] == TARGET.denominator, 'target denominator mismatch')
    for rel, digest in manifest['sha256'].items():
        path = (ROOT / rel).resolve()
        require(path.is_relative_to(ROOT), 'manifest path escapes package')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                f'hash mismatch: {rel}')
    parameters = json.loads((ROOT / 'parameters.json').read_text())
    require(parameters['segments'] == 160 and parameters['denominator'] == 10**10,
            'common parameter format mismatch')
    intervals = [{'source': f'proposals/t{2*j:06d}_w2.in', 'N': 10240,
                  'lo': 2*j, 'width': 2} for j in range(5120)]
    require({p.name for p in (ROOT / 'proposals').glob('*.in')} ==
            {Path(row['source']).name for row in intervals}, 'proposal coverage mismatch')
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
    for j, row in enumerate(intervals):
        require(row['lo'] == 2*j and row['lo']//64 == (row['lo']+1)//64,
                f'interval coverage mismatch: {j}')
    print('AUDIT PASS: all inputs, hashes, plane budgets, and 5,120 intervals.', flush=True)
    return intervals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'runs/recheck')
    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error('--jobs must be positive')
    intervals = audit()
    if args.audit_only:
        return
    output_dir = args.output_dir.resolve()
    require(not output_dir.exists() or not any(output_dir.iterdir()),
            'output directory must be empty')
    ids = list(range(len(intervals)))
    with tempfile.TemporaryDirectory(prefix='ps160-verify-') as work:
        work = Path(work)
        exe = work / 'verify'
        subprocess.run(['g++', '-O3', '-std=c++17',
                        str(ROOT / 'sources/strict160_verify.cpp'), '-o', str(exe)], check=True)
        def replay_shard(pair):
            k, selected = pair
            tasks, output = work / f'tasks{k}', work / f'output{k}.jsonl'
            tasks.write_text('\n'.join(f"{intervals[j]['source']} 10240 {intervals[j]['lo']} 2"
                                       for j in selected) + '\n')
            with output.open('w') as stream:
                subprocess.run([str(exe), str(tasks), str(TARGET.numerator),
                                str(TARGET.denominator), 'RECTANGLE_V1'],
                               cwd=ROOT, stdout=stream, check=True)
            got = [json.loads(line) for line in output.read_text().splitlines()]
            require(len(got) == len(selected), f'record count mismatch: shard {k}')
            for j, row in zip(selected, got):
                require(all(row[key] == intervals[j][key] for key in ('source','N','lo','width')),
                        f'computed interval assignment mismatch: {j}')
                require(row['status'] == 'PASS' and int(row['denominator']) == 2**56,
                        f'arithmetic verification failed: {j}')
                require(Fraction(int(row['lower_numerator']), int(row['denominator'])) >= TARGET,
                        f'computed lower bound below target: {j}')
            print(f'REPLAY PASS: shard {k}, {len(got)} intervals.', flush=True)
            return list(zip(selected, got))
        jobs = min(args.jobs, len(ids))
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            shards = list(pool.map(replay_shard, [(j, ids[j::jobs]) for j in range(jobs)]))
    rows = sorted(item for shard in shards for item in shard)
    require([j for j, _ in rows] == ids, 'computed arithmetic coverage is incomplete')
    records = [row for _, row in rows]
    minimum = min(Fraction(int(row['lower_numerator']), int(row['denominator'])) for row in records)
    require(minimum == EXPECTED_MINIMUM, 'computed minimum differs from the fixed certificate')
    report = {'status': 'VERIFIED', 'target_numerator': TARGET.numerator,
              'target_denominator': TARGET.denominator, 'target_met': minimum >= TARGET,
              'exact_target_margin': str(minimum-TARGET),
              'complete_input_and_coverage_audit': True,
              'all_5120_arithmetic_checks_replayed': True,
              'arithmetic_checks_replayed': len(records),
              'minimum_numerator': str(minimum.numerator),
              'minimum_denominator': str(minimum.denominator),
              'minimum_decimal': float(minimum)}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / 'records.json').write_text(json.dumps(records, indent=2) + '\n')
    (output_dir / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
