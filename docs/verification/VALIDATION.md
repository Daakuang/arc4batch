# Validation on 1 October 2026

The repository was checked from a separate copy containing only the files selected for upload.

- Eight numerical tests passed.
- All 24 retained scientific source files matched the new provenance record.
- The four NMPC engines passed an AST comparison with the evaluated source after module/path renaming. The serialized model hashes remain unchanged.
- The conditional tuning comparison checked 72,000 scalar steps and the independently specified zero-request actuator response.
- The reduced reference problems were solved again on 100- and 400-interval grids and reintegrated independently. The projected PI example and its derivative/refinement checks passed. The 400-interval reference conversion was 0.9223015870; the PI conversion was 0.9198908007.
- All 310 files listed in the archived-data manifest matched their checksums. Completion and recipe-conformity counts were recomputed from all 60 main trajectories.
- The summary, main trajectory plots, utility plots, and two LaTeX drawings were regenerated successfully. The main trajectory figure was visually inspected.

The full industrial simulation campaign was not rerun in this source cleanup. No controller parameters, process equations, stored trajectories, or scientific outcomes were changed. Platform-specific numerical differences remain possible.

Raw trajectories are supplied separately from the source package. Their checksum and file manifest support verified import and reanalysis.

The renamed MATLAB controller completed a two-step simulation with finite states; its overwrite guard rejected a repeated run. The renamed Python release reproduced the same complete 60-run summary, reduced demo, and figures from a clean copy.
