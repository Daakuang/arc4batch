"""Design-data parameter estimates and conditional prediction errors."""
import json
from artifact_paths import output_dir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adaptive_model import OUT,ROOT

def main():
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(3,2,figsize=(6.5,6.2),sharex=True,layout='constrained')
    metrics={}
    for j,scenario in enumerate(['N','PM_plus']):
        raw=json.loads((OUT/'adaptive_design'/f'observer_{scenario}_seed090.json').read_text())
        rows=raw['rows'];t=np.array([r['t'] for r in rows])/3600
        theta=np.array([r['adaptive']['theta'] for r in rows]);truth=[1,1] if j==0 else [1.2,.8]
        for k,(label,color,ls) in enumerate([(r'$\theta_H$','#1864A0','-'),(r'$\theta_{UA}$','#CE691C','--')]):
            axes[0,j].plot(t,theta[:,k],label=label,color=color,ls=ls,lw=1.3)
            axes[0,j].axhline(truth[k],color=color,ls=':',lw=1)
        for name,color,ls in [('nominal','#CE691C','--'),('adaptive','#218459','-')]:
            forecast=[r for r in rows if name+'_forecast_error' in r]
            tf=np.array([r['t'] for r in forecast])/3600
            e=np.array([r[name+'_forecast_error'] for r in forecast])
            for k in range(2):axes[k+1,j].plot(tf,e[:,k],color=color,ls=ls,lw=1,label={'nominal':'Nominal EKF','adaptive':'Augmented EKF'}[name])
        axes[0,j].set_title('N' if j==0 else r'PM$_+$')
        axes[2,j].set_xlabel('Time (h)')
        for k in range(3):
            axes[k,j].grid(alpha=.2);axes[k,j].text(.02,.95,f'({chr(97+k*2+j)})',transform=axes[k,j].transAxes,va='top')
            if k:axes[k,j].axhline(0,c='#555555',lw=.6)
        metrics[scenario]=raw['metrics']
    axes[0,0].set_ylabel('Parameter ratio')
    axes[1,0].set_ylabel('30 s reactor error (K)')
    axes[2,0].set_ylabel('30 s jacket error (K)')
    axes[0,0].legend(loc='lower right',ncol=2,fontsize=8)
    axes[1,0].legend(loc='lower right',fontsize=8)
    target=output_dir()/'parameter_validation'
    fig.savefig(target.with_suffix('.pdf'));fig.savefig(target.with_suffix('.png'),dpi=180);plt.close(fig)
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
