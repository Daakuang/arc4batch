"""Respect BOTH original benchmark thermal constraints with a derived selector."""
from artifact_paths import output_dir
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'results/jpc_revision_20260907/reduced'
import reduced_model as m
C=m.DH/(m.rho*m.cp);Qprime=2*m.U*(m.T-m.Tc_min)/m.sigma;kappa=.2

def law(y):
    xa,V,r=y;NA=m.Na0*(1-xa);NB=max(0,m.cb_feed*(V-m.Va0)-m.Na0*xa)
    rate=m.k_T*NA*NB/V;cb=NB/V
    e=m.Tmax-m.T-C*min(NA,NB)/V-m.e_sp
    uPI=np.clip(m.Kp*e+r,0,m.u_max)
    h=m.U*(2*V/m.sigma+np.pi*m.sigma**2)*(m.T-m.Tc_min)-m.DH*rate
    a=-m.DH*m.k_T*rate*(NA+NB)/V
    b=m.DH*m.k_T*NA/V*(m.cb_feed-cb)-Qprime
    uQ=np.clip((-a+kappa*h)/b,0,m.u_max) if b>0 else m.u_max
    u=min(uPI,uQ) if V<m.Vmax else 0.
    ri=m.Ki*e
    if (r>=m.r_plus and ri>0) or (r<=-m.r_minus and ri<0):ri=0.
    return u,ri,e,h,uPI,uQ,a,b,rate

def rhs(t,y):
    u,ri,*rest=law(y)
    return [rest[-1]/m.Na0,u,ri]

sol=solve_ivp(rhs,[0,m.tf],[0,m.Va0,0],rtol=1e-9,atol=[1e-11,1e-14,1e-15],dense_output=True,max_step=5)
assert sol.success
t=np.arange(0,m.tf+1);xa,V,r=sol.sol(t)
values=np.array([law(y) for y in np.array([xa,V,r]).T]);u=values[:,0];e=values[:,2];h=values[:,3]
NA=m.Na0*(1-xa);NB=m.cb_feed*(V-m.Va0)-m.Na0*xa
valid=(NA>=NB)&(NB>0)&(V<m.Vmax-1e-10)&(h>=-1e-6)
start=np.flatnonzero(valid&(e<=m.eta_e))[0]
end=np.flatnonzero((~valid)&(np.arange(len(t))>start))[0]
assert np.min(h)>=-1e-6,np.min(h)
assert np.max(m.Tmax-m.e_sp-e)<=m.Tmax+1e-6
assert V.max()<=m.Vmax+1e-10
assert e[start:end].min()>=-m.eta_s-1e-6 and e[start:end].max()<=m.eta_e+1e-6
indices=np.arange(start,end,17);dt=.01
reaction=m.k_T*NA[indices]*NB[indices]/V[indices]
newNA=NA[indices]-dt*reaction
newNB=NB[indices]+dt*(m.cb_feed*u[indices]-reaction)
newV=V[indices]+dt*u[indices]
eplus=m.Tmax-m.T-C*newNB/newV-m.e_sp
expected=C*m.k_T*NA[indices]*NB[indices]/V[indices]**2-C*(m.cb_feed-NB[indices]/V[indices])/V[indices]*u[indices]
derivative_error=float(np.max(np.abs((eplus-e[indices])/dt-expected)))
assert derivative_error<1e-6
refined=solve_ivp(rhs,[0,m.tf],[0,m.Va0,0],rtol=1e-10,atol=[1e-12,1e-15,1e-16],max_step=2)
assert refined.success
refinement_error=float(abs(refined.y[0,-1]-xa[-1]));assert refinement_error<1e-6
summary=json.loads((O/'summary.json').read_text())
cbe=(m.Tmax-m.T-m.e_sp-m.eta_e)/C
selector_margin=C*m.k_T*cbe**2
summary.update(ARC_conversion=float(xa[-1]),peak_Tcf_K=float((m.Tmax-m.e_sp-e).max()),
    controller_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    derivative_check_error=derivative_error,refined_conversion=float(refined.y[0,-1]),conversion_refinement_error=refinement_error,
    max_cooling_residual_W=float(-h.min()),regional_entry_h=t[start]/3600,regional_exit_h=t[end]/3600,
    e_min=float(e[start:end].min()),e_max=float(e[start:end].max()),kappa_Q_per_s=kappa,
    selector_upper_inward_lower_K_s=float(selector_margin),
    admissible_additive_residual_K_s=float(min(summary['safety_inward_lower_K_s'],selector_margin)),
    selector_activity_percent=float(100*np.mean((values[:,5]<values[:,4]-1e-14)&(V<m.Vmax-1e-10))),
    conversion_gap_percentage_points=float(100*(summary['references'][-1]['transcribed_conversion']-xa[-1])))
(O/'summary.json').write_text(json.dumps(summary,indent=2))
np.savez_compressed(O/'arc.npz',t=t,xa=xa,V=V,r=r,e=e,u=u,Tcf=m.Tmax-m.e_sp-e,
                    h_Q=h,u_PI=values[:,4],u_Q=values[:,5])
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':10,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'pdf.fonttype':42})
ref=np.load(O/'ocp_400.npz');fig,ax=plt.subplots(3,1,figsize=(5.4,6.2),layout='constrained')
ax[0].step(ref['t'][:-1]/3600,ref['u']*3.6e9,where='post',label='OCP reference (400 intervals)',color='#1864A0')
ax[0].plot(t/3600,u*3.6e9,'--',label='Projected PI + heat selector',color='#CE691C')
ax[0].set(ylabel='Feed (mL/h)',title='(a) Feed comparison');ax[0].legend()
ax[1].plot(t[start:end]/3600,e[start:end],color='#1864A0')
ax[1].axhline(-m.eta_s,c='black',ls='--');ax[1].axhline(m.eta_e,c='black',ls='--')
ax[1].set(ylabel='Error (K)',title='(b) First-order regional band')
ax[2].plot(t/3600,h,color='#218459');ax[2].axhline(0,c='black',ls='--')
ax[2].set(ylabel='Heat margin (W)',xlabel='Time (h)',title='(c) Instantaneous cooling capacity')
for a in ax:a.set_xlim(0,30);a.grid(alpha=.2)
fig.savefig(output_dir()/'reduced_benchmark_column.pdf')
fig.savefig(O/'reduced_benchmark.png',dpi=160)
print(json.dumps(summary,indent=2),flush=True)
