"""Paired NMPC evaluation aligned with the primary temperature-quality requirement."""
from pathlib import Path
import argparse,json,time,hashlib,sys
import numpy as np
import casadi as ca
from scipy.io import savemat
from parameter_model import ParameterModel,ParameterEKF,ROOT,OUT,MODELS,PLANT_OBS,MEAS_STD,SCALE_U,LOW_U,INIT_U,TAU,RATE,arr
from quality_controller import QualityController as ParameterController
from quality_model import QualityModel as ParameterModel

def run(a):
    controller='A_NMPC' if a.adaptive else 'NMPC'
    base=OUT/'quality_priority'/('evaluation' if a.evaluation else 'design')
    base.mkdir(parents=True,exist_ok=True)
    folder=base/f'{a.scenario}_{controller}_seed{a.seed:03d}'
    if (folder/'config.json').exists():raise RuntimeError('Refusing to overwrite a run: '+str(folder))
    folder.mkdir(exist_ok=True)
    sources=[Path(__file__),Path(__file__).with_name('parameter_model.py'),Path(__file__).with_name('quality_controller.py'),Path(__file__).with_name('quality_model.py')]
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    config=dict(vars(a),controller=controller,status='running',source_hashes=source_hashes,
                source_sha256=hashlib.sha256(json.dumps(source_hashes,sort_keys=True).encode()).hexdigest(),
                T_ref=351.,T_constraint=351.7,T_recipe_lower=350.3,T_recipe_hold_start_s=3300,T_trip=373.15,
                objective_temperature_weight=a.temperature_weight,objective_common_scale=.001,
                priority='two-sided temperature recipe first during t>=3300 s; strong hold-temperature tracking; pressure target secondary; original startup stage weights retained',
                nominal_parameters=[327240.,21496.,363.,283.],input_lags=TAU.tolist(),input_rates=RATE.tolist(),
                estimator_period_s=1,controller_period_s=30,prediction_horizon_s=a.horizon*30,
                parameter_initial=[1.,1.],parameter_bounds=[.5,1.5],initial_parameter_std=[.2,.2],
                parameter_random_walk_std_per_s=a.parameter_std,
                actuator_prediction='initial version: reduced process predictor; common physical lag/rate-limited plant and actual applied inputs in EKF',
                operating_constraints='hard-first two-sided hold temperature recipe and pressure upper bound; strongly penalized recovery; all true departures retained',
                solver_latency='offline synchronous solve; measured time does not advance plant clock')
    (folder/'config.json').write_text(json.dumps(config,indent=2))
    zn=np.random.Generator(np.random.PCG64(a.seed)).standard_normal((5,30001))
    noise_path=OUT/f'noise_{a.seed:03d}.mat'
    if not noise_path.exists():savemat(noise_path,dict(paired_noise=zn))
    pdata=json.loads((MODELS/f'plant_{int(a.scenario=="F")}.json').read_text())
    x=arr(pdata['x0']);d=arr(pdata['d0']);u=INIT_U.copy()
    if a.scenario in ['PM_plus','F']:d[:2]*=[1.2,.8]
    if a.scenario=='PM_minus':d[:2]*=[.8,1.2]
    plant=ca.Function.load(str(MODELS/f'plant_flow_{int(a.scenario=="F")}.casadi'))
    model=ParameterModel(a.temperature_weight);observer=ParameterEKF(model,a.adaptive,a.parameter_std)
    control=ParameterController(model,a.horizon)
    timer=time.perf_counter();states=[x.copy()];inputs=[];commands=[];estimates=[];obslog=[];solverlog=[]
    reason='time_limit';command=INIT_U.copy()
    for t in range(a.max_steps):
        y=x[PLANT_OBS]+MEAS_STD*zn[:3,t]
        obs=observer.update(y,x[[15,16]])
        if t%30==0:
            command,info=control.control(observer.z,t)
            info['t']=t;solverlog.append(info);obslog.append(dict(t=t,**obs))
            estimates.append(np.r_[t,observer.z[:12]*model.scale,observer.z[12:]])
            print(json.dumps(dict(T=float(x[11]),dose=float(x[15]),theta=observer.z[12:].tolist(),**info)),flush=True)
            (folder/'solver_progress.json').write_text(json.dumps(solverlog,indent=2))
            if t and t%600==0:
                np.savez_compressed(folder/'checkpoint.npz',t=np.arange(len(states)),x=np.array(states),
                    u=np.array(inputs),u_cmd=np.array(commands),estimate=np.array(estimates))
        u=np.clip(u+np.clip((np.clip(command,LOW_U,SCALE_U)-u)/TAU,-RATE,RATE),LOW_U,SCALE_U)
        u[:2]=np.minimum(u[:2],np.maximum(0,np.array([3250.,.002])-x[[15,16]]))
        commands.append(command.copy());inputs.append(u.copy())
        x=arr(plant(x0=x,p=np.r_[u,d])['xf']);states.append(x.copy());observer.predict(u)
        if x[11]>=373.15:reason='temperature_trip';break
        if x[15]>=3250-1e-3 and x[13]<=1e6:reason='completed';break
        if not np.isfinite(x).all():reason='nonfinite_plant';break
    config.update(status='finished',termination=reason,steps=len(inputs),actual_plant_parameters=d.tolist(),
                  wall_seconds=time.perf_counter()-timer)
    np.savez_compressed(folder/'trajectory.npz',t=np.arange(len(states)),x=np.array(states),u=np.array(inputs),
                        u_cmd=np.array(commands),estimate=np.array(estimates))
    (folder/'solver_log.json').write_text(json.dumps(solverlog,indent=2))
    (folder/'observer_log.json').write_text(json.dumps(obslog,indent=2))
    (folder/'config.json').write_text(json.dumps(config,indent=2))
    print(json.dumps(config),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scenario',required=True,choices=['N','PM_minus','PM_plus','F'])
    p.add_argument('--seed',type=int,default=91);p.add_argument('--adaptive',action='store_true')
    p.add_argument('--evaluation',action='store_true');p.add_argument('--max-steps',type=int,default=30000)
    p.add_argument('--temperature-weight',type=float,default=3e6/360);
    p.add_argument('--horizon',type=int,default=60);p.add_argument('--parameter-std',type=float,default=1e-4)
    run(p.parse_args())
