"""Run the examples and reproduction commands from this repository root."""
from pathlib import Path
import argparse
import runpy
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'src' / 'arc4batch'


def run_script(name, arguments=()):
    """Keep the numerical engines' original import and command-line interfaces."""
    sys.path.insert(0, str(SOURCE))
    sys.argv = [name + '.py', *arguments]
    runpy.run_path(str(SOURCE / (name + '.py')), run_name='__main__')


def import_results(archive):
    """Import only verified result files from the matching companion ZIP."""
    import hashlib
    import json

    archive = Path(archive).resolve()
    asset = json.loads((ROOT / 'docs/verification/release-assets.json').read_text())['results']
    if hashlib.sha256(archive.read_bytes()).hexdigest() != asset['sha256']:
        raise ValueError('This archive does not match the recorded result asset.')
    with zipfile.ZipFile(archive) as source:
        planned = []
        for entry in source.infolist():
            parts = Path(entry.filename).parts
            if 'results' not in parts or entry.is_dir():
                continue
            relative = Path(*parts[parts.index('results'):])
            target = (ROOT / relative).resolve()
            if not target.is_relative_to(ROOT / 'results'):
                raise ValueError('Unexpected archive path: ' + entry.filename)
            content = source.read(entry)
            if target.exists() and target.read_bytes() != content:
                raise FileExistsError('Refusing to replace a different result: ' + str(relative))
            planned.append((entry, target))
        for entry, target in planned:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read(entry))
    print('Imported the matching archived results without replacing different runs.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['check', 'demo', 'import-results', 'summarize', 'plot', 'noise', 'nmpc', 'suite', 'tuning'])
    args, remaining = parser.parse_known_args()
    if args.command == 'check':
        sys.path.insert(0, str(ROOT / 'tools'))
        sys.argv = ['verify_release.py', *remaining]
        runpy.run_path(str(ROOT / 'tools/verify_release.py'), run_name='__main__')
    elif args.command == 'demo':
        run_script('run_reduced_reference')
        run_script('reduced_with_heat_selector')
    elif args.command == 'import-results':
        if len(remaining) != 1:
            parser.error('import-results needs the path to the matching results ZIP')
        import_results(remaining[0])
    elif args.command == 'summarize':
        run_script('analyze_results', ['--require-complete', *remaining])
    elif args.command == 'plot':
        run_script('plot_revision_figures', remaining or ['trajectories'])
    elif args.command == 'tuning':
        sys.path.insert(0, str(SOURCE))
        import json
        from vpc_tuning import screen
        print(json.dumps([screen(gain, integral_time) for gain, integral_time in [(6,40),(8,80),(8,40),(10,40)]], indent=2))
    else:
        target = {'noise':'prepare_noise', 'nmpc':'run_quality_nmpc', 'suite':'run_quality_suite'}[args.command]
        run_script(target, remaining)


if __name__ == '__main__':
    main()
