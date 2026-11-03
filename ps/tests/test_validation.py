"""Regression tests for validation under ordinary and optimized Python."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='ps validation ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for source in ROOT.iterdir():
            dest = self.root / source.name
            if source.name in {'replay.py', 'finalize_manifest.py',
                               'certificate.json', 'RESULTS.json'}:
                shutil.copy2(source, dest)
            else:
                dest.symlink_to(source, target_is_directory=source.is_dir())

    def run_modes(self, script, args, error):
        for label, flags, optimize in (
                ('normal', [], None), ('-O', ['-O'], None),
                ('PYTHONOPTIMIZE', [], '1')):
            with self.subTest(mode=label):
                env = os.environ.copy()
                env.pop('PYTHONOPTIMIZE', None)
                if optimize is not None:
                    env['PYTHONOPTIMIZE'] = optimize
                result = subprocess.run(
                    [sys.executable, *flags, str(self.root / script), *args],
                    cwd=self.root, env=env, capture_output=True, text=True,
                    check=False)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(error, result.stderr)
                self.assertNotIn('AUDIT PASS', result.stdout)

    def test_wrong_target_rejected(self):
        path = self.root / 'certificate.json'
        manifest = json.loads(path.read_text())
        manifest['target_numerator'] = 1
        path.write_text(json.dumps(manifest))
        self.run_modes('replay.py', ['--audit-only'], 'target numerator mismatch')

    def test_bad_checksum_rejected(self):
        path = self.root / 'certificate.json'
        manifest = json.loads(path.read_text())
        manifest['sha256']['parameters.json'] = '0' * 64
        path.write_text(json.dumps(manifest))
        self.run_modes('replay.py', ['--audit-only'], 'hash mismatch: parameters.json')

    def test_nonpositive_jobs_rejected(self):
        self.run_modes('replay.py', ['--audit-only', '--jobs', '0'],
                       '--jobs must be positive')

    def test_manifest_refresh_rejects_invalid_result(self):
        path = self.root / 'RESULTS.json'
        result = json.loads(path.read_text())
        result['target_met'] = False
        path.write_text(json.dumps(result))
        before = (self.root / 'certificate.json').read_bytes()
        self.run_modes('finalize_manifest.py', [], 'unexpected certified result')
        self.assertEqual((self.root / 'certificate.json').read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
