"""Reproduce the 100/400-interval reference and independently replay its feed."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.integrate import solve_ivp
import reduced_model as m
ROOT=Path(__file__).resolve().parents[2]
dest=ROOT/'results/jpc_revision_20260907/reduced';dest.mkdir(parents=True,exist_ok=True)

def main():
    C=m.DH/(m.rho*m.cp)
    cs=(m.Tmax-m.T-m.e_sp+m.eta_s)/C
    ce=(m.Tmax-m.T-m.e_sp-m.eta_e)/C
    us=float(np.clip(-m.Kp*m.eta_s+m.r_plus,0,m.u_max))
    ue=float(np.clip(m.Kp*m.eta_e-m.r_minus,0,m.u_max))
    qmin=C*(m.cb_feed-ce)/m.Vmax
    required=m.k_T*m.Na0*ce/(m.cb_feed-ce)
    lower=C*m.k_T*cs**2;upper=qmin*(ue-required)
    assert us==0 and lower>0 and upper>0
    summary=dict(order=1,Kp=m.Kp,Ki=m.Ki,r_plus=m.r_plus,r_minus=m.r_minus,
                 alpha_s_max=us,alpha_e_min=ue,c_B_s=cs,c_B_e=ce,
                 required_feed_upper=required,q_e_lower=qmin,
                 safety_inward_lower_K_s=lower,economic_inward_lower_K_s=upper,
                 admissible_additive_residual_K_s=min(lower,upper),
                 model_sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest())
    refs=[]
    for N in [100,400]:
        tt,aa,vv,uu=m.solve_ocp(N)
        x=np.array([0.,m.Va0]);maxcf=m.T;maxq=-np.inf
        for j in range(N):
            def rhs(_,y):return [m.rate(y[0],y[1])*y[1]/m.Na0,uu[j]]
            check=solve_ivp(rhs,[tt[j],tt[j+1]],x,rtol=1e-10,atol=1e-13,max_step=20)
            assert check.success
            x=check.y[:,-1]
            for a,v in check.y.T:
                maxcf=max(maxcf,m.Tcf(a,v))
                maxq=max(maxq,m.DH*m.rate(a,v)*v-m.U*(2*v/m.sigma+np.pi*m.sigma**2)*(m.T-m.Tc_min))
        refs.append(dict(N=N,numerical_guard_K=.001,numerical_guard_W=.001,
            transcribed_conversion=float(aa[-1]),independent_conversion=float(x[0]),
            peak_Tcf_K=float(maxcf),max_cooling_residual_W=float(maxq),
            final_dynamics_error=float(np.max(np.abs(x-np.r_[aa[-1],vv[-1]])))))
        np.savez_compressed(dest/f'ocp_{N}.npz',t=tt,xa=aa,V=vv,u=uu)
        print(json.dumps(refs[-1]),flush=True)
    summary['references']=refs
    (dest/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
