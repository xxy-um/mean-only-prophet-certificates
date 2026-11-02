#!/usr/bin/env python3
"""Rebuild the 5120 final certificate proposals from frozen rational seeds.

This reproduces the deterministic final assembly, not the historical search
which originally selected the seed threshold and plane lists.  Python 3.9+,
standard library only.
"""

from pathlib import Path
from fractions import Fraction
from collections import defaultdict, Counter
import argparse
import bisect
import hashlib
import json

ROOT = Path(__file__).resolve().parent
MASS_DENOMINATOR = 2**32
N = 10240
WIDTH = 2


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rebuild(output_dir, compare_dir=None):
    output_dir = Path(output_dir).resolve()
    if compare_dir is not None:
        compare_dir = Path(compare_dir).resolve()
        if output_dir == compare_dir:
            raise ValueError('output directory must differ from the comparison directory')
    output_dir.mkdir(parents=True, exist_ok=True)

    seed_hashes = json.loads((ROOT / 'SEED_SHA256SUMS.json').read_text())
    for name, expected in seed_hashes.items():
        if digest(ROOT / name) != expected:
            raise ValueError('seed checksum mismatch: ' + name)
    expected_hashes = json.loads((ROOT / 'FINAL_PROPOSAL_SHA256SUMS.json').read_text())
    seeds = sorted(json.loads((ROOT / 'plane_choice_seeds.json').read_text()),
                   key=lambda row: Fraction(row['lo']))
    starts = [Fraction(row['lo']) for row in seeds]
    parameters = json.loads((ROOT / 'parameters.json').read_text())
    common = (parameters['b_numerators'] + parameters['m_numerators']
              + parameters['kappa_numerators'] + [0] * 160)
    if parameters['segments'] != 160 or parameters['denominator'] != 10**10:
        raise ValueError('unexpected common-parameter format')
    cache = {}

    def read_seed(name):
        if name in cache:
            return cache[name]
        lines = (ROOT / name).read_text().splitlines()
        if lines[0] != 'MEAN_ONLY_AFFINE_V1':
            raise ValueError('seed proposal format: ' + name)
        coarse, segments = map(int, lines[1].split())
        if segments != 160 or list(map(int, lines[2].split())) != common:
            raise ValueError('seed common-parameter mismatch: ' + name)
        tail = list(map(int, lines[3].split()))
        if len(tail) != coarse or len(lines) - 5 != int(lines[4]):
            raise ValueError('seed proposal dimensions: ' + name)
        planes = []
        for line in lines[5:]:
            fields = line.split()
            planes.append((int(fields[0]), int(fields[1]), tuple(fields[2:])))
        result = coarse, lines[2], tail, planes
        cache[name] = result
        return result

    found = {}
    byte_comparisons = 0
    for lo in range(0, N, WIDTH):
        midpoint = Fraction(2 * lo + WIDTH, 2 * N)
        index = bisect.bisect_right(starts, midpoint) - 1
        if index < 0:
            raise ValueError('seed interval does not cover midpoint')
        row = seeds[index]
        if not Fraction(row['lo']) <= midpoint <= Fraction(row['hi']):
            raise ValueError('seed interval does not cover midpoint')
        left = Fraction(row['anchor_left'])
        right = Fraction(row['anchor_right'])
        alpha = (midpoint - left) / (right - left)
        if not 0 <= alpha <= 1:
            raise ValueError('mixture coefficient out of range')

        c1, params1, tail, p1 = read_seed(row['source_left'])
        c2, params2, _, p2 = read_seed(row['source_right'])
        if params1 != params2:
            raise ValueError('mixture uses different common parameters')
        coarse = max(c1, c2)
        if coarse % c1 or coarse % c2:
            raise ValueError('seed owner grids are not nested')
        masses = defaultdict(int)
        for count, planes, weight in ((c1, p1, 1 - alpha), (c2, p2, alpha)):
            refinement = coarse // count
            for owner, mass, rest in planes:
                kind, cell = divmod(owner, count)
                weighted_mass = mass * weight.numerator // weight.denominator
                if not weighted_mass:
                    continue
                for child in range(refinement):
                    masses[(kind * coarse + cell * refinement + child, rest)] += weighted_mass
        totals = Counter()
        for (owner, _), mass in masses.items():
            totals[owner] += mass
        if max(totals.values(), default=0) > MASS_DENOMINATOR:
            raise ValueError('owner mixture mass exceeds one')

        lines = ['MEAN_ONLY_AFFINE_V1', f'{coarse} 160', params1,
                 ' '.join(str(value) for value in tail for _ in range(coarse // c1)),
                 str(len(masses))]
        for (owner, rest), mass in sorted(masses.items()):
            lines.append(' '.join([str(owner), str(mass), *rest]))
        raw = ('\n'.join(lines) + '\n').encode('utf-8')
        name = f't{lo:06d}_w{WIDTH}.in'
        actual = hashlib.sha256(raw).hexdigest()
        if expected_hashes.get(name) != actual:
            raise ValueError('rebuilt proposal checksum mismatch: ' + name)
        if compare_dir is not None:
            if (compare_dir / name).read_bytes() != raw:
                raise ValueError('rebuilt proposal byte mismatch: ' + name)
            byte_comparisons += 1
        (output_dir / name).write_bytes(raw)
        found[name] = actual

    if found != expected_hashes:
        raise ValueError('final proposal coverage mismatch')
    manifest_bytes = json.dumps(found, sort_keys=True, separators=(',', ':')).encode()
    report = {
        'status': 'PASS',
        'generation_scope': 'deterministic final assembly from frozen rational seeds',
        'segments': 160,
        'N_p': N,
        'N_t': N // WIDTH,
        'seed_intervals': len(seeds),
        'seed_proposals': len(cache),
        'final_proposals': len(found),
        'expected_proposal_hashes_matched': len(found),
        'byte_comparisons_passed': byte_comparisons,
        'proposal_manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
    }
    (output_dir / 'REBUILD_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT.parent / 'regenerated_proposals')
    parser.add_argument('--compare-dir', type=Path,
                        help='optionally compare every generated file byte-for-byte with this directory')
    args = parser.parse_args()
    print(json.dumps(rebuild(args.output_dir, args.compare_dir), indent=2))
