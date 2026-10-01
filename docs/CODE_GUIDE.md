# Reading the implementation

The process model, controllers, simulation drivers, and reporting are separate. The root `arc4batch.py` only selects a task; it does not contain control equations.

## Reduced reactor

Read `src/arc4batch/reduced_model.py` first. It defines physical constants, the reaction rate, the temperature after cooling failure, and the guarded reference optimal-control problem. `solve_reduced_reference.py` solves the reference and checks it with a separate numerical integrator. `reduced_arc.py` then computes the projected PI command, limits it by instantaneous cooling capacity, integrates the closed loop, and checks the error-coordinate derivative.

The PI error is a thermal coordinate; it is not the measured reactor temperature error. Its sign and allowable band are defined by the reduced model. Both the reference feed and the regulated feed have the same physical feed and volume limits.

## Industrial ARC

The physical reactor equations are in `matlab/ptfe.m`. Start reading the controller at `matlab/arc_controller.m`: measured temperature drives the thermal controller; its demand passes through the split-range utility map; virtual demand before limits drives the economic PI; that command changes initiator feed and the pressure setpoint. Dose limits and the actuator dynamics determine the actual inputs applied to the plant.

`simulate_arc.m` runs one scenario. `run_arc_suite.m` supplies the paired scenario/seed combinations. `run_arc_additional.m` and `run_tuning_study.m` run the separate stress and tuning studies. `arc_controller_stress.m` contains the prescribed actuator-lag stress variant. `export_models.m` exports the CasADi functions used by Python.

## NMPC and state estimation

1. `adaptive_model.py` loads the industrial model, sets the physical input scaling, and defines the EKF. The adaptive EKF estimates reaction-heat and heat-transfer multipliers.
2. `nmpc_model.py` defines the common temperature-priority objective and operating constraints.
3. `nmpc_controller.py` solves the prediction problem, checks the accepted solution, and applies the declared recovery or fallback action when required.
4. `simulate_nmpc.py` applies measurements, estimation, optimization, actuator dynamics, and plant integration in that order. Estimation updates every second and NMPC every 30 seconds.

These four engines retain the evaluated numerical expressions, with renamed modules and result paths. `run_nmpc_suite.py` verifies the release hashes before running paired batches. `prepare_noise.py` generates the common seeded measurement noise. `run_nmpc_stress.py` changes only the declared noise and actuator specification.

## Results and figures

`analyze_results.py` reads the true trajectories and evaluates completion, temperature-recipe conformity, operating limits, and tracking windows. A temperature trip is not a completed batch. A completed batch that violates the recipe is retained in the record but excluded from the primary conforming-batch time comparison.

`plot_results.py` draws the study trajectories and utility signals. `plot_parameter_validation.py` draws the estimator checks. All generated figures go to `outputs/`. The two physical/recipe drawings use the small LaTeX sources in `docs/figure_sources/` and need a LaTeX installation only when regenerated.

## What is deliberately excluded

Reviewer replies, marked manuscripts, draft table generators, exploratory runs, dependency binaries, temporary logs, and the working repository's history are not part of this code repository. The optional data archive retains both successful and failed scientific runs; those outcomes are evidence rather than temporary clutter.
