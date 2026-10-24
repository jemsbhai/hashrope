"""Compare the historical and released exact engine after rustfmt only."""
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    historical = ROOT / 'benchmarks/exact_optimization_round4/optimized.zip'
    baseline = HERE / 'baseline.zip'
    report = {
        'historical_archive_sha256': digest(historical),
        'baseline_archive_sha256': digest(baseline),
        'tool': subprocess.check_output(['rustfmt', '--version'], text=True).strip(),
        'scope': 'comparison section RangeStart through separator before Free-function wrappers; formatting only, no assembly claim',
        'script_sha256': digest(Path(__file__)),
    }
    sections = {}
    with tempfile.TemporaryDirectory(prefix='hashrope-comparison-source-') as directory:
        for label, archive in [('historical', historical), ('baseline', baseline)]:
            with zipfile.ZipFile(archive) as zipped:
                source = zipped.read('hashrope/packages/rust/src/rope.rs').decode().replace('\r\n', '\n')
            start = source.index('#[derive(Clone, Copy)]\nstruct RangeStart')
            end = source.index('Free-function wrappers', start)
            target = Path(directory) / (label + '.rs')
            target.write_text(source[start:end], encoding='utf-8')
            subprocess.run(['rustfmt', '--edition', '2021', str(target)], check=True)
            sections[label] = target.read_bytes().replace(b'\r\n', b'\n')
            report[label + '_formatted_comparison_sha256'] = hashlib.sha256(sections[label]).hexdigest()
    report['comparison_sections_identical_after_rustfmt'] = sections['historical'] == sections['baseline']
    assert report['comparison_sections_identical_after_rustfmt']
    # Print rather than overwrite the original separately recorded audit.
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
