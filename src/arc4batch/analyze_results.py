"""True-state metrics with explicit windows, completion and failure accounting."""
from pathlib import Path
from artifact_paths import output_dir
import json,csv,re,argparse
import numpy as np
from scipy.io import loadmat
import casadi as ca
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results/jpc_revision_20260907'
FINAL_SEEDS=range(200,205)
CONTROLLERS=['ARC','NMPC','A_NMPC']

def run_folder(scenario,controller,seed):
    base=OUT if controller.startswith('ARC') else OUT/'quality_priority/evaluation'
    return base/f'{scenario}_{controller}_seed{seed:03d}'

def stress_folder(scenario,controller):
    base=OUT/'stress' if controller.startswith('ARC') else OUT/'quality_priority/stress'
    return base/f'{scenario}_{controller}_seed050'

def metric(t,x,u,*,scenario,controller,seed,solver=None,declared_status=None):
    t=np.asarray(t);x=np.asarray(x);u=np.asarray(u)
    assert len(t)==len(x) and len(u)==len(t)-1 and np.all(np.diff(t)==1)
    assert np.isfinite(x).all() and np.isfinite(u).all()
    temp=x[:,11];pressure=x[:,13];dose=x[:,15]
    trip=bool(np.max(temp)>=373.15)
    complete=bool(dose[-1]>=3250-1e-3 and pressure[-1]<=1e6 and not trip)
    hold=t>=3300
    feeding=hold&(dose<3250-1e-3)
    def errors(values,mask):
        if not np.any(mask):return [None,None]
        e=np.abs(values[mask]);return [float(np.mean(e)),float(np.max(e))]
    et=errors(temp-351.,hold);ep=errors((pressure-1.5e6)/1000,feeding)
    result=dict(scenario=scenario,controller=controller,seed=int(seed),duration_h=float(t[-1]/3600),
       completed=complete,termination='completed' if complete else ('temperature_trip' if trip else 'incomplete'),
       T_MAE_K=et[0],T_absmax_K=et[1],P_MAE_kPa=ep[0],P_absmax_kPa=ep[1],
       T_peak_K=float(np.max(temp)),T_excess_K=float(max(0,np.max(temp)-351.7)),
       T_exceedance_percent=float(100*np.mean(temp[1:]>351.7)),
       T_exceedance_seconds=int(np.count_nonzero(temp[1:]>351.7)),
       P_peak_MPa=float(np.max(pressure)/1e6),P_exceedance_seconds=int(np.count_nonzero(pressure[1:]>1.6e6)),
       final_dose_A_kg=float(dose[-1]),final_dose_B_kg=float(x[-1,16]),final_pressure_MPa=float(pressure[-1]/1e6),
       samples=len(t),hold_samples=int(np.count_nonzero(hold)),pressure_samples=int(np.count_nonzero(feeding)))
    result['T_recipe_band_exceedance_percent']=float(100*np.mean(np.abs(temp[hold]-351)>.7)) if np.any(hold) else None
    result['T_below_recipe_percent']=float(100*np.mean(temp[hold]<350.3)) if np.any(hold) else None
    result['T_recipe_departure_samples']=int(np.sum(np.abs(temp[hold]-351)>.7))
    result['temperature_recipe_conforming']=bool(np.any(hold) and result['T_recipe_departure_samples']==0)
    result['operating_conforming']=bool(result['T_exceedance_seconds']==0 and result['P_exceedance_seconds']==0)
    result['conforming_completed']=bool(complete and result['temperature_recipe_conforming'] and result['operating_conforming'])
    # Independent moment calculation, using the final full-model state.
    den0=x[-1,0]+x[-1,3];den1=x[-1,1]+x[-1,4]
    mn=100*den1/den0 if den0>0 else np.nan
    mw=100*(x[-1,2]+x[-1,5])/den1 if den1>0 else np.nan
    result.update(Mn_g_mol=float(mn),Mw_g_mol=float(mw),PDI=float(mw/mn),
                  polymer_mass_kg=float(.1*den1))
    if solver is not None:
        times=np.array([v['seconds'] for v in solver]);ok=np.array([v['accepted'] for v in solver])
        result.update(solver_calls=len(solver),solver_rejected=int(np.sum(~ok)),
           solver_not_converged=sum(v['status'] not in ['Solve_Succeeded','Solved_To_Acceptable_Level'] for v in solver),
           solve_median_s=float(np.median(times)),solve_p95_s=float(np.quantile(times,.95)),
           solve_max_s=float(np.max(times)),solves_over_30s=int(np.sum(times>30)))
    if declared_status is not None:
        result['declared_termination']=declared_status
        assert (declared_status=='completed')==complete, (controller,scenario,seed,declared_status,result)
    # Operating slew applies before the ideal metered dose-stop interlock.
    assert np.max(u[:,0])<=.694+1e-9 and np.max(u[:,1])<=4e-7+1e-15
    assert np.min(u[:,:2])>=-1e-12 and np.max(np.abs(u[:,2]))<=1+1e-9
    assert np.max(dose)<=3250+1e-5 and np.max(x[:,16])<=.002+1e-10
    return result

def load_run(folder):
    match=re.match(r'(N|PM_minus|PM_plus|F)_(.+)_seed(\d+)$',folder.name)
    if not match:raise ValueError(folder.name)
    scenario,controller,seed=match.groups();seed=int(seed)
    if controller.startswith('ARC'):
        data=loadmat(folder/'trajectory.mat',simplify_cells=True);s=data['SimData']
        if 'u_cmd_actual' not in s:raise RuntimeError('superseded ARC signal version '+folder.name)
        x=np.column_stack((s['initial_state'],s['X'])).T
        u=s['u_impl'].T;t=np.arange(len(x))
        row=metric(t,x,u,scenario=scenario,controller=controller,seed=seed)
        row['wall_seconds']=float(data['elapsed'])
        act=(t[:-1]>=3300)&(x[:-1,15]<3250-1e-3)
        row['virtual_min']=float(np.min(s['alpha'][act])) if np.any(act) else None
        row['Psp_min_MPa']=float(np.min(s['Psp'][act])/1e6) if np.any(act) else None
        row['Psp_relax_percent']=float(100*np.mean(s['Psp'][act]<1.5e6-1)) if np.any(act) else None
        row['FB_active_percent']=float(100*np.mean(u[act,1]>1e-12)) if np.any(act) else None
        row['cooling_overload_percent']=float(100*np.mean(s['alpha'][act]<=-1)) if np.any(act) else None
        row['Psp_all_min_MPa']=float(np.min(s['Psp'])/1e6)
        row['solve_median_s']=float(np.median(s['solve_times']))
        row['solve_p95_s']=float(np.quantile(s['solve_times'],.95))
        row['solve_max_s']=float(np.max(s['solve_times']))
        signals=dict(Psp=s['Psp'],virtual=s['alpha'],integral=s['ui_VPC'],economic=s['val_FB'],commands=s['u_cmd_actual'].T)
    else:
        cfg=json.loads((folder/'config.json').read_text())
        if cfg['status']!='finished':raise RuntimeError('unfinished '+folder.name)
        data=np.load(folder/'trajectory.npz');t=data['t'];x=data['x'];u=data['u']
        log=json.loads((folder/'solver_log.json').read_text())
        row=metric(t,x,u,scenario=scenario,controller=controller,seed=seed,solver=log,declared_status=cfg['termination'])
        row['wall_seconds']=cfg['wall_seconds'];row['source_sha256']=cfg.get('source_sha256')
        signals=dict(estimate=data['estimate'],commands=data['u_cmd'])
    return row,dict(t=t,x=x,u=u,**signals)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--require-complete',action='store_true');a=parser.parse_args()
    records=[];pending=[]
    for seed in FINAL_SEEDS:
        for scenario in ['N','PM_minus','PM_plus','F']:
            for controller in CONTROLLERS:
                name=f'{scenario}_{controller}_seed{seed:03d}';folder=run_folder(scenario,controller,seed)
                try: row,_=load_run(folder);records.append(row)
                except (FileNotFoundError,RuntimeError):pending.append(name)
    fields=sorted(set().union(*(r.keys() for r in records))) if records else []
    with (output_dir()/'run_metrics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fields);writer.writeheader();writer.writerows(records)
    groups=[]
    for scenario in ['N','PM_minus','PM_plus','F']:
        for controller in CONTROLLERS:
            rows=[r for r in records if r['scenario']==scenario and r['controller']==controller]
            done=[r for r in rows if r['completed']]
            conforming=[r for r in rows if r['conforming_completed']]
            group=dict(scenario=scenario,controller=controller,n=len(rows),completed=len(done),
                       conforming_completed=len(conforming),
                       seeds=[r['seed'] for r in rows],incomplete_seeds=[r['seed'] for r in rows if not r['completed']])
            times=[r['duration_h'] for r in conforming]
            for suffix,fn in [('mean',np.mean),('min',np.min),('max',np.max)]:
                group['conforming_duration_h_'+suffix]=float(fn(times)) if times else None
            for key in ['duration_h','Mn_g_mol','Mw_g_mol','PDI','polymer_mass_kg']:
                values=[r[key] for r in done if r.get(key) is not None]
                group[key+'_mean']=float(np.mean(values)) if values else None
                group[key+'_min']=float(np.min(values)) if values else None
                group[key+'_max']=float(np.max(values)) if values else None
            for key in ['T_MAE_K','P_MAE_kPa']:
                values=[r[key] for r in rows if r.get(key) is not None]
                group[key+'_mean']=float(np.mean(values)) if values else None
                group[key+'_min']=float(np.min(values)) if values else None
                group[key+'_max']=float(np.max(values)) if values else None
            group['temperature_window_runs']=sum(r['hold_samples']>0 for r in rows)
            group['pressure_window_runs']=sum(r['pressure_samples']>0 for r in rows)
            for key in ['T_peak_K','T_absmax_K','P_absmax_kPa','T_excess_K','T_exceedance_percent','T_recipe_band_exceedance_percent','T_below_recipe_percent','P_peak_MPa']:
                values=[r[key] for r in rows if r.get(key) is not None]
                group[key+'_max']=float(np.max(values)) if values else None
            group['violating_runs']=sum(r['T_exceedance_seconds']>0 for r in rows)
            group['recipe_departing_runs']=sum(r['T_recipe_departure_samples']>0 for r in rows)
            group['pressure_violating_runs']=sum(r['P_exceedance_seconds']>0 for r in rows)
            group['solver_rejected_total']=sum(r.get('solver_rejected',0) for r in rows)
            groups.append(group)
    summary=dict(expected_runs=60,available_runs=len(records),pending=pending,groups=groups)
    (output_dir()/'results_summary.json').write_text(json.dumps(summary,indent=2))
    if a.require_complete:
        assert not pending,pending
        protocol=json.loads((OUT/'quality_priority/frozen_protocol.json').read_text())
        expected=protocol['source_sha256']
        assert all(r['source_sha256']==expected for r in records if r['controller']!='ARC')
        for r in records:
            if r['controller']=='ARC':continue
            cfg=json.loads((run_folder(r['scenario'],r['controller'],r['seed'])/'config.json').read_text())
            assert cfg['temperature_weight']==protocol['settings']['temperature_weight']
            assert cfg['horizon']*30==protocol['settings']['horizon_s']
            assert cfg['parameter_std']==protocol['settings']['parameter_random_walk_std']
    print(json.dumps(dict(available=len(records),pending=len(pending),completed=sum(r['completed'] for r in records))))

if __name__=='__main__':main()
