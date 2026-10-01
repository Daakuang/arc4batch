"""Physical parameter augmentation shared by the nominal and adaptive NMPC baselines."""
from pathlib import Path
import json
import casadi as ca
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/study'
MODELS=OUT/'models'
SCALE_U=np.array([.694,4e-7,1.])
LOW_U=np.array([0.,0.,-1.])
INIT_U=np.array([.694,0.,0.])
TAU=np.array([5.,5.,3.]);RATE=.05*SCALE_U
MAP=np.array([0,6,7,8,9,11,12,13,14,15,16,18])
OBS=np.array([5,8,7]);PLANT_OBS=np.array([11,14,13])
MEAS_STD=np.array([.1/3,.1/3,5000/3])

def arr(v):return np.array(v,dtype=float).reshape(-1)

class ParameterModel:
    def __init__(self):
        p=json.loads((MODELS/'reduced.json').read_text())
        self.d=arr(p['d0']);self.scale=arr(p['x_norm']);self.scale[10]=.002
        self.initial=arr(p['x0'])/self.scale
        self.lo=np.array([0,0,0,0,0,320.15,0,1.5e5,283,0,0,2000])/self.scale
        self.hi=np.array([np.inf]*5+[373.15,10000,2e6,373,3250,.002,6000])/self.scale
        red=ca.Function.load(str(MODELS/'reduced_f.casadi'))
        z=ca.SX.sym('z',12);u=ca.SX.sym('u',3);th=ca.SX.sym('theta',2)
        dp=ca.vertcat(self.d[0]*th[0],self.d[1]*th[1],self.d[2:])
        dx,cost=red(z*self.scale,u*SCALE_U,dp)
        hot=.5*(u[2]+ca.sqrt(u[2]**2+1e-6));cold=.5*(-u[2]+ca.sqrt(u[2]**2+1e-6))
        x=z*self.scale
        dx[8]=(10*4200*cold*(self.d[3]-x[8])+4200*hot*(self.d[2]-x[8])
               -dp[1]*(x[8]-x[5]))/8400000
        self.f=ca.Function('physical_parameter_rhs',[z,u,th],[dx/self.scale,cost])
        za=ca.MX.sym('za',14);ua=ca.MX.sym('ua',3)
        rr=self.f(za[:12],ua,za[12:])[0]
        flow=ca.integrator('parameter_observer_flow','idas',
            dict(x=za,p=ua,ode=ca.vertcat(rr,ca.MX.zeros(2))),0,1,
            dict(abstol=1e-9,reltol=1e-7,max_num_steps=10000))
        zp=flow(x0=za,p=ua)['xf']
        self.flow=ca.Function('parameter_state_flow',[za,ua],[zp])
        self.step=ca.Function('parameter_observer_step',[za,ua],[zp,ca.jacobian(zp,za)])

class ParameterEKF:
    def __init__(self,model,adaptive,parameter_std=1e-4,noise_multiplier=1.):
        self.model=model;self.adaptive=adaptive
        self.z=np.r_[model.initial,1.,1.]
        scale=np.r_[model.scale,1.,1.]
        sd=np.array([1e-9,1e-5,.02,1e-8,1e-4,.002,.001,10.,.002,0.,0.,.001,
                     parameter_std if adaptive else 0.,parameter_std if adaptive else 0.])
        init=np.array([1e-7,1e-4,1.,1e-6,.05,.03,.01,1000.,.03,0.,0.,.1,
                       .2 if adaptive else 0.,.2 if adaptive else 0.])
        self.Q=np.diag((sd/scale)**2);self.P=np.diag((init/scale)**2)
        self.R=np.diag((MEAS_STD*noise_multiplier/model.scale[OBS])**2)
        self.H=np.eye(14)[OBS]

    def update(self,y,doses):
        y=np.asarray(y)/self.model.scale[OBS]
        innovation=y-self.H@self.z
        S=self.H@self.P@self.H.T+self.R
        K=np.linalg.solve(S,(self.P@self.H.T).T).T
        self.z+=K@innovation
        IK=np.eye(14)-K@self.H
        self.P=IK@self.P@IK.T+K@self.R@K.T
        self.z[:12]=np.maximum(self.z[:12],self.model.lo)
        self.z[[9,10]]=np.asarray(doses)/self.model.scale[[9,10]]
        self.P[[9,10],:]=0;self.P[:,[9,10]]=0
        before=self.z[12:].copy()
        self.z[12:]=np.clip(self.z[12:],.5,1.5) if self.adaptive else 1.
        return dict(innovation=(innovation*self.model.scale[OBS]).tolist(),
                    NIS=float(innovation@np.linalg.solve(S,innovation)),
                    theta=self.z[12:].tolist(),parameter_clipped=bool(np.any(before!=self.z[12:])),
                    parameter_std=np.sqrt(np.maximum(0,np.diag(self.P)[12:])).tolist())

    def predict(self,u):
        zp,A=self.model.step(self.z,np.asarray(u)/SCALE_U)
        self.z=arr(zp);A=np.array(A)
        self.P=A@self.P@A.T+self.Q;self.P=(self.P+self.P.T)/2
        if not self.adaptive:self.z[12:]=1.
