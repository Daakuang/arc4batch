# Reproducing the study

Run commands from the repository root. Use the recorded Python dependencies and MATLAB/CasADi versions in `verification/environment-original.txt` for the closest numerical reproduction. Solver trajectories and timings can differ across platforms.

## Quick example and checks

```sh
python arc4batch.py check
python arc4batch.py demo
python arc4batch.py tuning
```

The reduced example requires no archived industrial results. It writes generated outputs to `outputs/` and its numeric records to `results/jpc_revision_20260907/reduced/`. Use a separate copy if you want to preserve a previously imported archive unchanged.

## Archived results

The matching data archive is `arc4batch-v1.0.0-results.zip`; its checksum is in `verification/release-assets.json`. The 1.0.1 source cleanup reuses this unchanged result asset.

```sh
python arc4batch.py import-results /path/to/arc4batch-v1.0.0-results.zip
python arc4batch.py check --results
python arc4batch.py summarize
python arc4batch.py plot
python arc4batch.py plot signals
```

`check --results` verifies every file in the data manifest and recomputes completion/conformity counts from all 60 main trajectories. `summarize` writes `outputs/run_metrics.csv` and `outputs/results_summary.json`. `plot` generates the main trajectory comparison; `plot signals` generates the physical utility and virtual-demand figure.

The recorded directory name `jpc_revision_20260907` is retained inside `results/` because archived configurations and hashes refer to it. Readers do not need to navigate that directory to use the root commands.

## Rerun an industrial case

Full industrial batches can take substantial time. Use a separate reproduction copy; runners refuse to overwrite completed or incompatible existing runs.

```sh
python arc4batch.py noise --seeds 200
python arc4batch.py nmpc --scenario N --seed 200 --evaluation
python arc4batch.py nmpc --scenario N --seed 200 --evaluation --adaptive
```

In MATLAB, set `root` to the repository root and add the `matlab` folder:

```matlab
root = pwd;
addpath(fullfile(root, 'matlab'));
run_arc_revision(root, 'N', 200, 30000, 'full_arc');
```

The MATLAB CasADi distribution must be available in `.runtime/casadi/matlab`, as expected by the runners. It is not bundled. Saved CasADi functions in the source repository support Python reproduction without another MATLAB export.

## Full main comparison

Five paired seeds, 200–204, four scenarios, and three controllers give 60 main runs. `N` is nominal; `PM_plus` and `PM_minus` have the declared opposite reaction-heat/heat-transfer mismatches; `F` adds the gel-effect fault to the adverse mismatch.

```sh
python arc4batch.py noise
python arc4batch.py suite --workers 2
```

```matlab
run_arc_suite(root, 200:204, false);
```

The suite preserves the frozen controller parameters, paired measurements, optimizer acceptance checks, and physical actuator/dose recurrences. Results computed under different source hashes are not combined. No simulations impose computation delay on the plant clock.

## Additional evidence

The supplied data archive contains the separate ablation, tuning, observer, and stress records. The main comparison is still 60 runs; additional design batches are not pooled into it.

```matlab
run_arc_suite(root, 0, true);
run_arc_additional(root);
run_industrial_tuning_bridge(root);
```

```sh
python src/arc4batch/run_quality_stress.py --scenario N
python src/arc4batch/run_quality_stress.py --scenario PM_plus
```

For the local economic-loop checks, `python arc4batch.py tuning` evaluates the four declared gain/integral-time pairs. Conditions and response deadlines are in `config/vpc_tuning.json`. The deadline for a pressure setpoint is not a guarantee for the actual pressure response.

## Figures and model export

```sh
python arc4batch.py plot diagrams
python src/arc4batch/plot_parameter_validation.py
```

Regenerating the two LaTeX drawings requires `latexmk` and the usual TikZ/PGFPlots packages. Other plots use Matplotlib. In MATLAB, `export_models(root)` regenerates the industrial CasADi models; doing so changes model provenance and should be done in a separate copy.

## Verification record

`verification/source-provenance.json` maps each retained source to its validated predecessor. Four frozen NMPC source files are byte-identical. Metric calculations and tuning functions retain their numerical expressions. The cleanup changes entry points, artifact destinations, figure labels, MATLAB search paths, and documentation.

Before making this repository public, attach the matching result archive and record the final validation result. The repository remains private until the author requests public release.
