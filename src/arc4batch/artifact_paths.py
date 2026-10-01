"""Create the public export's generated-artifact directory when needed."""
from pathlib import Path


def output_dir():
    target = Path(__file__).resolve().parents[2] / 'outputs'
    target.mkdir(parents=True, exist_ok=True)
    return target
