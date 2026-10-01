"""Check source provenance, frozen models, and optionally all archived data."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/arc4batch'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', action='store_true', help='Also check the complete imported results archive.')
    args = parser.parse_args()
    provenance = json.loads((ROOT/'docs/verification/source-provenance.json').read_text())
    for entry in provenance['files']:
        assert digest(ROOT/entry['path']) == entry['sha256'], entry['path']
    from run_nmpc_suite import verify_source
    from adaptive_model import ParameterModel
    from vpc_tuning import audit_envelopes, screen
    protocol = verify_source()
    ParameterModel()
    assert screen(8,40)['passes'] and not screen(8,80)['passes']
    envelopes = audit_envelopes()
    tests = unittest.defaultTestLoader.discover(str(ROOT/'verification'))
    result = unittest.TextTestRunner(verbosity=1).run(tests)
    assert result.wasSuccessful()
    report = {'source_files_verified':len(provenance['files']), 'numerical_tests':result.testsRun,
              'release_source_and_evaluated_models':'verified', 'conditional_tuning_checks':envelopes}
    if args.results:
        manifest = json.loads((ROOT/'docs/verification/results-manifest.json').read_text())
        for entry in manifest:
            path = ROOT/entry['path']
            assert path.is_file(), 'Import the matching results ZIP first: '+entry['path']
            assert digest(path) == entry['sha256'], entry['path']
        from analyze_results import OUT, load_run, run_folder
        summary = json.loads((OUT/'results_summary.json').read_text())
        runs = 0
        for group in summary['groups']:
            rows = [load_run(run_folder(group['scenario'],group['controller'],seed))[0] for seed in group['seeds']]
            assert sum(row['completed'] for row in rows) == group['completed']
            assert sum(row['conforming_completed'] for row in rows) == group['conforming_completed']
            runs += len(rows)
        assert runs == 60
        report.update(archived_files_verified=len(manifest), main_runs_recomputed=runs)
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
