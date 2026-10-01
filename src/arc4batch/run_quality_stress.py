"""Synthetic stress with frozen quality-priority settings and unchanged EKF covariance."""
import argparse,hashlib,json
from pathlib import Path
import run_quality_nmpc as runner
from run_quality_suite import verify_source


if __name__=='__main__':
    protocol=verify_source()
    p=argparse.ArgumentParser();p.add_argument('--scenario',required=True,choices=['N','PM_plus']);a=p.parse_args()
    a.seed=50;a.adaptive=True;a.evaluation=True;a.max_steps=30000;a.horizon=60;a.parameter_std=1e-4
    a.temperature_weight=protocol['settings']['temperature_weight']
    a.noise_multiplier=3;a.lag_multiplier=2
    a.wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    root=(runner.OUT/'quality_priority/stress').resolve();root.mkdir(parents=True,exist_ok=True)
    name=f'{a.scenario}_A_NMPC_seed050';destination=(root/name).resolve()
    assert destination.is_relative_to(root) and not destination.exists(),destination
    runner.OUT=root;runner.MEAS_STD=runner.MEAS_STD*3;runner.TAU=runner.TAU*2
    source=(root/'quality_priority/evaluation'/name).resolve()
    assert source.is_relative_to(root)
    if not (source/'config.json').exists():runner.run(a)
    config=json.loads((source/'config.json').read_text())
    assert config['status']=='finished' and config['source_sha256']==protocol['source_sha256']
    assert config['wrapper_sha256']==a.wrapper_sha256
    # Rehome this one verified completed case within its named stress directory;
    # preserve all files and refuse overwrite. No recursive deletion is used.
    source.rename(destination)
    print(json.dumps(dict(completed_stress_directory=str(destination))))
