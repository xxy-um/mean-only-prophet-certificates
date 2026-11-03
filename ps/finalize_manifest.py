#!/usr/bin/env python3
"""Refresh integrity metadata after adding completed audit records."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json

ROOT = Path(__file__).resolve().parent
for required in ('README.md', 'RELEASE_NOTES.md',
                 'Verification.pdf'):
    if not (ROOT / required).is_file():
        raise FileNotFoundError('missing release file: ' + required)
result = json.loads((ROOT / 'RESULTS.json').read_text())
if not (result['target_met'] and result['minimum_numerator'] == '49839735470020217'):
    raise ValueError('unexpected certified result')
target = Fraction(result['target_numerator'], result['target_denominator'])
minimum = Fraction(int(result['minimum_numerator']), int(result['minimum_denominator']))
if target != Fraction(1729, 2500) or minimum < target:
    raise ValueError('unexpected target or minimum below target')
if Fraction(result['exact_target_margin']) != minimum-target:
    raise ValueError('incorrect exact target margin')
manifest = dict(result)
manifest['sha256'] = {
    p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(ROOT.rglob('*'))
    if p.is_file() and p.name != 'certificate.json'
    and not {'.git', '__pycache__', 'regenerated_proposals'}.intersection(p.relative_to(ROOT).parts)
    and not p.name.endswith('.pyc')
    and p.name not in {'Verification.aux', 'Verification.log',
                       'Verification.out', 'Verification.synctex.gz'}
}
(ROOT / 'certificate.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
print(f"Manifest records {len(manifest['sha256'])} files.")
