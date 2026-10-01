"""Economic VPC withdrawal and recovery bounds for the declared operating region.

screen(K, Ti) checks the response deadlines in config/vpc_tuning.json.
These local action bounds require the sustained demand conditions stated there;
they do not certify the complete reactor's thermal invariance.
"""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
P = json.loads((ROOT / 'config/vpc_tuning.json').read_text())
C = P['controller']
REQ = P['chosen_engineering_requirements']

def actuator(b,cmd):
    return np.clip(b+np.clip((cmd-b)/C['initiator_tau_s'],
             -C['initiator_slew_kg_s2'],C['initiator_slew_kg_s2']),0,C['initiator_max_kg_s'])

def commands(v):
    return (np.clip(v*C['initiator_gain_kg_s'],0,C['initiator_max_kg_s']),
            np.clip(C['pressure_nominal_Pa']+C['pressure_gain_Pa']*np.minimum(v,0),
                    C['pressure_min_Pa'],C['pressure_nominal_Pa']))

def envelopes(K,Ti,mode,steps=300):
    eta=.15 if mode=='overload' else .25
    n=np.arange(1,steps+1)
    if mode=='overload':
        r=np.maximum(C['integral_min'],C['integral_max']-n*K/Ti*eta)
        v=-K*eta+r; b=C['initiator_max_kg_s']
    else:
        r=np.minimum(C['integral_max'],C['integral_min']+n*K/Ti*eta)
        v=K*eta+r; b=0.
    v=np.clip(v,*C['economic_output_bounds']);cmd,psp=commands(v)
    actual=[]
    for c in cmd:
        b=float(actuator(b,c));actual.append(b)
    return dict(r=r,v=v,cmd=cmd,psp=psp,actual=np.array(actual))

def first_time(condition):
    ids=np.flatnonzero(condition)
    return int(ids[0]+1) if len(ids) else None

def screen(K,Ti):
    s=envelopes(K,Ti,'overload');e=envelopes(K,Ti,'spare')
    ptime=first_time(s['psp']<=C['pressure_nominal_Pa']-REQ['pressure_target_reduction_Pa']+1e-8)
    ptime=None if ptime is None else ptime+C['pressure_execution_delay_samples']
    btime=first_time(e['actual']>=.5*C['initiator_max_kg_s']-1e-18)
    immediate=bool(K*.15>=C['integral_max']-1e-12)
    return dict(K=K,Ti_s=Ti,immediate_zero_command_sufficient=immediate,
                pressure_target_deadline_bound_s=ptime,half_feed_recovery_bound_s=btime,
                initiator_one_percent_bound_s=first_time(s['actual']<=.01*C['initiator_max_kg_s']+1e-18),
                residual_initiator_dose_bound_kg=float(np.sum(s['actual'])),
                passes=bool(immediate and ptime is not None and ptime<=30 and btime is not None and btime<=60))

def audit_envelopes():
    # Independent scalar loop with random demand and starting integral/actuator.
    rng=np.random.default_rng(7301);checks=0
    for K in P['analytic_grid']['K']:
        for Ti in P['analytic_grid']['Ti_s']:
            for mode in ['overload','spare']:
                bound=envelopes(K,Ti,mode,100)
                for trial in range(40):
                    r=float(rng.uniform(-1.5,1));b=float(rng.uniform(0,4e-7))
                    for n in range(100):
                        err=(-.15-rng.uniform(0,2)) if mode=='overload' else (.25+rng.uniform(0,1))
                        r=min(1,max(-1.5,r+K/Ti*err))
                        v=min(10,max(-40,K*err+r))
                        cmd=min(4e-7,max(0,1e-7*v))
                        # Separately written recurrence, including slew and saturation.
                        diff=(cmd-b)/5
                        b=min(4e-7,max(0,b+min(2e-8,max(-2e-8,diff))))
                        if mode=='overload': assert b<=bound['actual'][n]+1e-18 and r<=bound['r'][n]+1e-12
                        else: assert b>=bound['actual'][n]-1e-18 and r>=bound['r'][n]-1e-12
                        checks+=1
    # Hand result: 15 rate-limited steps then geometric decay with ratio 0.8.
    b=4e-7;values=[]
    for i in range(500):b=float(actuator(b,0));values.append(b)
    assert abs(values[14]-1e-7)<1e-20
    assert first_time(np.array(values)<=4e-9)==30
    assert abs(sum(values)-4e-6)<1e-18
    return dict(random_comparison_samples=checks,random_seed=7301,
                zero_command_1percent_s=30,zero_command_remaining_dose_kg=4e-6)
