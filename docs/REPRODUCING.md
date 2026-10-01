# Reproducing the study

Run commands from the repository root. Use the recorded Python dependencies and MATLAB/CasADi versions in `verification/environment-original.txt` for the closest numerical reproduction. Solver trajectories and timings can differ across platforms.

## Quick example and checks

```sh
python arc4batch.py check
python arc4batch.py demo
python arc4batch.py tuning
```

The reduced example requires no archived industrial results. It writes generated outputs to `outputs/` and its numeric records to `results/study/reduced/`. Use a separate copy if you want to preserve a previously imported archive unchanged.

## Archived results

The matching data archive is `arc4batch-v1.0.2-results.zip`; its checksum is in `verification/release-assets.json`. This archive uses the same descriptive paths as the source release. Raw trajectories and run metadata retain their evaluated contents.

```sh
python arc4batch.py import-results /path/to/arc4batch-v1.0.2-results.zip
python arc4batch.py check --results
python arc4batch.py summarize
python arc4batch.py plot
python arc4batch.py plot signals
```

`check --results` verifies every file in the data manifest and recomputes completion/conformity counts from all 60 main trajectories. `summarize` writes `outputs/run_metrics.csv` and `outputs/results_summary.json`. `plot` generates the main trajectory comparison; `plot signals` generates the physical utility and virtual-demand figure.

Recorded and generated study results are grouped under `results/study/`; figures and summaries go to `outputs/`.

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
simulate_arc(root, 'N', 200, 30000, 'full_arc');
```

The MATLAB CasADi distribution must be available in `.runtime/casadi/matlab`, as expected by the runners. It is not bundled. Saved CasADi functions in the source repository support Python reproduction without another MATLAB export.

## Full main comparison

Five paired seeds, 200â€“204, four scenarios, and three controllers give 60 main runs. `N` is nominal; `PM_plus` and `PM_minus` have the declared opposite reaction-heat/heat-transfer mismatches; `F` adds the gel-effect fault to the adverse mismatch.

```sh
python arc4batch.py noise
python arc4batch.py suite --workers 2
```

```matlab
run_arc_suite(root, 200:204, false);
```

The suite preserves the evaluated controller parameters, paired measurements, optimizer acceptance checks, and physical actuator/dose recurrences. It recognizes the evaluated source and this release's module/path renaming; unrecognized source versions are rejected. No simulations impose computation delay on the plant clock.

## Additional evidence

The supplied data archive contains the separate ablation, tuning, observer, and stress records. The main comparison is still 60 runs; additional design batches are not pooled into it.

```matlab
run_arc_suite(root, 0, true);
run_arc_additional(root);
run_tuning_study(root);
```

```sh
python src/arc4batch/run_nmpc_stress.py --scenario N
python src/arc4batch/run_nmpc_stress.py --scenario PM_plus
```

For the local economic-loop checks, `python arc4batch.py tuning` evaluates the four declared gain/integral-time pairs. Conditions and response deadlines are in `config/vpc_tuning.json`. The deadline for a pressure setpoint is not a guarantee for the actual pressure response.

## Figures and model export

```sh
python arc4batch.py plot diagrams
python src/arc4batch/plot_parameter_validation.py
```

Regenerating the two LaTeX drawings requires `latexmk` and the usual TikZ/PGFPlots packages. Other plots use Matplotlib. In MATLAB, `export_models(root)` regenerates the industrial CasADi models; doing so changes model provenance and should be done in a separate copy.

## Verification record

`verification/source-provenance.json` maps each retained source to its validated predecessor. An AST comparison verifies that the four NMPC engines differ only in module names and artifact paths. Metric calculations and tuning functions retain their numerical expressions. The cleanup changes module names, entry points, artifact destinations, figure labels, MATLAB search paths, and documentation. Both evaluation and release source identities are recorded.

The source package supports the reduced example and new industrial simulations. The companion archive additionally supports immediate reanalysis of the evaluated study.
