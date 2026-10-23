# arc4batch

Advanced regulatory control (ARC) for cooling-limited exothermic semi-batch reactors, with a reduced example and an industrial polymerization model.

The temperature loop requests cooling. An economic valve position controller uses the available cooling margin to adjust initiator feed and the pressure setpoint. The code also implements nominal and parameter-adaptive output-feedback NMPC for the paper's comparison.

**Published paper:** [A theory-guided advanced regulatory control synthesis for cooling-limited exothermic semi-batch reactors](https://doi.org/10.1016/j.jprocont.2026.103872), Chenchen Zhou and José Matias. *Journal of Process Control* **168** (2026), Article **103872**.

**DOI:** [10.1016/j.jprocont.2026.103872](https://doi.org/10.1016/j.jprocont.2026.103872) · **Preprint:** [arXiv:2606.18799](https://arxiv.org/abs/2606.18799)

**Version:** 1.0.2. Research code accompanying the published paper.

## How to cite

If you use the ARC method, code or simulation results in your research, please cite the published paper:

> Chenchen Zhou and José Matias. A theory-guided advanced regulatory control synthesis for cooling-limited exothermic semi-batch reactors. *Journal of Process Control* **168** (2026), Article **103872**. https://doi.org/10.1016/j.jprocont.2026.103872

```bibtex
@article{ZhouMatias2026ARC,
  author  = {Zhou, Chenchen and Matias, Jos{\'e}},
  title   = {A theory-guided advanced regulatory control synthesis for cooling-limited exothermic semi-batch reactors},
  journal = {Journal of Process Control},
  volume  = {168},
  pages   = {103872},
  year    = {2026},
  doi     = {10.1016/j.jprocont.2026.103872},
  url     = {https://doi.org/10.1016/j.jprocont.2026.103872}
}
```

The DOI above identifies the journal article. For reproducibility, also record the software release or commit used; the current release is [v1.0.2](https://github.com/Daakuang/arc4batch/releases/tag/v1.0.2). GitHub's **Cite this repository** menu uses the published paper specified in [CITATION.cff](CITATION.cff).

## Start with the reduced example

Use Python 3.12 and install the recorded dependencies in a virtual environment:

```sh
python -m pip install -r requirements.txt
python arc4batch.py check
python arc4batch.py demo
```

The small `verification/` directory holds the reproducibility checks, rather than development utilities. `check` runs eight numerical tests, verifies the release source and evaluated model hashes, and checks the conditional economic-loop response bounds. `demo` solves the 100- and 400-interval reduced reference problems, independently integrates the optimized feed, and simulates the projected PI controller with its cooling-capacity selector. It requires no MATLAB installation. The figure is written to `outputs/reduced_benchmark_column.pdf` and the numerical record to `results/study/reduced/`.

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

## Related papers and original model

The industrial reactor model is shared with these earlier studies:

- **Model and fault-monitoring benchmark:** Simin Li, Shuang-hua Yang, Yi Cao, Xiaoping Jiang and Chenchen Zhou, *A benchmark of industrial polymerization process for thermal runaway process monitoring*, Process Safety and Environmental Protection 193 (2025), 353–363. [DOI: 10.1016/j.psep.2024.11.057](https://doi.org/10.1016/j.psep.2024.11.057).
- **Real-time NMPC:** Chenchen Zhou, Zuzhen Ji and Jose Matias, *Real-time nonlinear model predictive control framework for event-triggered switching in industrial batch polymerization process*, Journal of Process Control 163 (2026), 103738. [DOI: 10.1016/j.jprocont.2026.103738](https://doi.org/10.1016/j.jprocont.2026.103738) · [Preprint](https://arxiv.org/abs/2606.14976).
- **Original public Simulink model:** Simin Li, [A Semi-batch Polymerization Process for Fault Simulation](https://www.mathworks.com/matlabcentral/fileexchange/169341-a-semi-batch-polymerization-process-for-fault-simulation), MATLAB Central File Exchange; first released in July 2024.

This repository accompanies the ARC paper. It includes the shared reactor equations and the nominal/adaptive NMPC implementations used for its comparison. It does not provide a complete reproduction package for the earlier real-time NMPC framework or the original fault-monitoring study. The File Exchange distribution is a separate Simulink implementation with its own license.
