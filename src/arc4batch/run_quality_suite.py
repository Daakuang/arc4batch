"""Execute the frozen temperature-priority paired evaluation."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,os,subprocess,sys,time
from parameter_model import ROOT,OUT


def verify_source():
    p=json.loads((OUT/'quality_priority/frozen_protocol.json').read_text())
    for name,digest in p['source_hashes'].items():
        assert hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()==digest,name
    for name,digest in p['model_hashes'].items():
        assert hashlib.sha256((OUT/'models'/name).read_bytes()).hexdigest()==digest,name
    return p


def run_case(task):
    scenario,seed,adaptive=task;protocol=verify_source()
    name=f'{scenario}_{"A_NMPC" if adaptive else "NMPC"}_seed{seed:03d}'
    cfg=OUT/'quality_priority/evaluation'/name/'config.json'
    if cfg.exists():
        c=json.loads(cfg.read_text())
        if (c.get('status')=='finished' and c['source_sha256']==protocol['source_sha256']
            and c['temperature_weight']==protocol['settings']['temperature_weight']):
            return dict(run=name,status='existing')
        raise RuntimeError('Existing unfinished or different-version run: '+name)
    command=[sys.executable,str(Path(__file__).with_name('run_quality_nmpc.py')),
             '--scenario',scenario,'--seed',str(seed),'--evaluation',
             '--temperature-weight',str(protocol['settings']['temperature_weight'])]
    if adaptive:command.append('--adaptive')
    env=os.environ.copy();env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
    start=time.perf_counter();logs=OUT/'quality_priority/logs';logs.mkdir(exist_ok=True)
    with (logs/f'{name}.log').open('w',encoding='utf-8') as f:
        result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
    return dict(run=name,status='finished' if result.returncode==0 else 'execution_error',
                returncode=result.returncode,seconds=time.perf_counter()-start)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=10);a=p.parse_args()
    protocol=verify_source();tasks=[(s,n,c=='A_NMPC') for n in protocol['seeds']
                                  for s in protocol['scenarios'] for c in ['NMPC','A_NMPC']]
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run_case,t) for t in tasks]):
            r=future.result();results.append(r);print(json.dumps(r),flush=True)
            (OUT/'quality_priority/suite_status.json').write_text(json.dumps(results,indent=2))
    assert all(r['status'] in ['finished','existing'] for r in results),results
