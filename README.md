# arc4batch

Advanced regulatory control (ARC) for cooling-limited exothermic semi-batch reactors, with a reduced example and an industrial polymerization model.

The temperature loop requests cooling. An economic valve position controller uses the available cooling margin to adjust initiator feed and the pressure setpoint. The code also implements nominal and parameter-adaptive output-feedback NMPC for the paper's comparison.

**Paper:** [A Theory-Guided Advanced Regulatory Control Synthesis for Cooling-Limited Exothermic Semi-Batch Reactors](https://arxiv.org/abs/2606.18799), Chenchen Zhou and Jose Matias.

**Version:** 1.0.2. Research code accompanying the linked preprint.

## Start with the reduced example

Use Python 3.12 and install the recorded dependencies in a virtual environment:

```sh
python -m pip install -r requirements.txt
python arc4batch.py check
python arc4batch.py demo
```

`check` runs eight numerical tests, verifies the release source and evaluated model hashes, and checks the conditional economic-loop response bounds. `demo` solves the 100- and 400-interval reduced reference problems, independently integrates the optimized feed, and simulates the projected PI controller with its cooling-capacity selector. It requires no MATLAB installation. The figure is written to `outputs/reduced_benchmark_column.pdf` and the numerical record to `results/study/reduced/`.

## Find the code you need

| Task | Start here |
|---|---|
| Understand the reduced reactor and PI example | `src/arc4batch/reduced_model.py` and `reduced_arc.py` |
| Read the industrial ARC controller | `matlab/arc_controller.m` |
| Read the industrial physical model | `matlab/ptfe.m` |
| Read the adaptive model and EKF | `src/arc4batch/adaptive_model.py` |
| Read the NMPC objective and constraints | `src/arc4batch/nmpc_model.py` and `nmpc_controller.py` |
| Understand completion and temperature metrics | `src/arc4batch/analyze_results.py` |
| Check VPC withdrawal and recovery bounds | `src/arc4batch/vpc_tuning.py` |
| Reproduce the full study | [docs/REPRODUCING.md](docs/REPRODUCING.md) |

The numerical engines use descriptive module names. Their numerical expressions match the evaluated version; imports and result paths have been updated. [docs/CODE_GUIDE.md](docs/CODE_GUIDE.md) explains their inputs, outputs, and execution order.

## Work with the archived study

The source repository includes serialized models and frozen configurations. The optional full trajectories are a separate, approximately 313 MB archive, `arc4batch-v1.0.2-results.zip`. The source repository contains no bulk trajectories. To use the separately supplied companion archive:

```sh
python arc4batch.py import-results /path/to/arc4batch-v1.0.2-results.zip
python arc4batch.py check --results
python arc4batch.py summarize
python arc4batch.py plot
```

Import verifies the archive checksum and refuses to replace different existing results. Summary and figure files go to `outputs/`; the original result files remain in their recorded structure. The main comparison contains 60 prescribed runs. Additional ablation, tuning, and stress runs remain separate.

MATLAB with CasADi is needed to rerun the industrial ARC simulations. Python and the saved CasADi model functions suffice for the NMPC comparison. The reproduction guide gives the commands and explains the archive layout.

## Scope and attribution

The tests and examples reproduce the stated models and scenarios. The local tuning bounds rely on the conditions in `config/vpc_tuning.json`; they do not certify the full industrial reactor. Batch-time comparisons include only completed runs that conform to the temperature recipe and operating limits. Recorded failures are retained.

Code and supplied simulated results use the [MIT license](LICENSE). Software author and maintainer: Chenchen Zhou. Cite [CITATION.cff](CITATION.cff) and the scientific paper when using the method or results. See [THIRD_PARTY.md](THIRD_PARTY.md) for dependency attribution.
