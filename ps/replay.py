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


def audit():
    manifest = json.loads((ROOT / 'certificate.json').read_text())
    assert manifest['target_numerator'] == TARGET.numerator
    assert manifest['target_denominator'] == TARGET.denominator
    for rel, digest in manifest['sha256'].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, ('hash mismatch', rel)
    parameters = json.loads((ROOT / 'parameters.json').read_text())
    intervals = json.loads((ROOT / 'intervals.json').read_text())
    records = json.loads((ROOT / 'records.json').read_text())
    assert len(intervals) == len(records) == 5120
    expected = None
    for source in sorted({x['source'] for x in intervals}):
        rows = (ROOT / source).read_text().splitlines()
        assert rows[0] == 'MEAN_ONLY_AFFINE_V1'
        coarse, k = map(int, rows[1].split())
        assert k == 160 and coarse > 0 and coarse % 160 == 0 and 10240 % coarse == 0
        params = tuple(map(int, rows[2].split()))
        assert len(params) == 640
        if expected is None:
            expected = params
        assert params == expected, ('different common parameters', source)
        b, m, kap, eta = (params[j:j+160] for j in range(0, 640, 160))
        assert all(v > 0 for v in b) and all(x < y for x, y in zip(b, b[1:]))
        assert sum(b) <= 160 * 10**10
        assert all(0 <= x <= y for x, y in zip(m, b))
        assert all(x >= 0 for x in kap) and all(x == 0 for x in eta)
        tail = list(map(int, rows[3].split()))
        assert len(tail) == coarse and tail == sorted(tail)
        assert all(0 <= x <= 2**20 for x in tail)
        count = int(rows[4])
        assert len(rows) == count + 5
        mass = [0] * (3 * coarse)
        for line in rows[5:]:
            v = line.split()
            owner, weight = int(v[0]), int(v[1])
            assert 0 <= owner < 3 * coarse and 0 <= weight <= 2**32
            mass[owner] += weight
            if v[2] == 'log':
                assert len(v) == 6 and owner < 2 * coarse and int(v[3]) == 1
            elif v[2] == 'C':
                assert len(v) == 4 and owner >= 2 * coarse and 0 <= int(v[3]) <= 10**10
            else:
                raise AssertionError(('non-article plane', v[2]))
        assert all(x <= 2**32 for x in mass)
    assert list(expected[:160]) == parameters['b_numerators']
    assert list(expected[160:320]) == parameters['m_numerators']
    assert list(expected[320:480]) == parameters['kappa_numerators']
    assert Fraction(sum(expected[:160]), 160 * 10**10) == Fraction(parameters['theta_0'])
    for j, (x, y) in enumerate(zip(intervals, records)):
        assert x['N'] == 10240 and x['lo'] == 2*j and x['width'] == 2
        assert x['lo'] // 64 == (x['lo'] + 1) // 64
        assert all(x[k] == y[k] for k in ('source', 'N', 'lo', 'width'))
        assert y['status'] == 'PASS'
        assert Fraction(int(y['lower_numerator']), int(y['denominator'])) >= TARGET
    minimum = min(Fraction(int(y['lower_numerator']), int(y['denominator'])) for y in records)
    assert minimum == Fraction(int(manifest['minimum_numerator']), int(manifest['minimum_denominator']))
    print(f'AUDIT PASS: target={TARGET} ({float(TARGET)}), all 5120 intervals, common parameters, plane budgets, hashes and coverage.', flush=True)
    return intervals, records, minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--critical-only', action='store_true')
    args = parser.parse_args()
    assert args.jobs > 0
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
            tasks.write_text('\n'.join(f"{ROOT / intervals[j]['source']} 10240 {intervals[j]['lo']} 2" for j in selected) + '\n')
            with output.open('w') as stream:
                subprocess.run([str(exe), str(tasks), str(TARGET.numerator), str(TARGET.denominator), 'RECTANGLE_V1'], stdout=stream, check=True)
            got = [json.loads(line) for line in output.read_text().splitlines()]
            assert len(got) == len(selected)
            for j, row in zip(selected, got):
                assert row['status'] == 'PASS'
                for key in KEYS:
                    assert row[key] == records[j][key], ('arithmetic replay mismatch', j, key)
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
