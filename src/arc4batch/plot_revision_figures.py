"""Scientific vector figures at the manuscript's actual 5.4-inch text width."""
import string
from artifact_paths import output_dir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from analyze_results import load_run,ROOT,OUT,run_folder
D=output_dir()
COLORS={'ARC':'#1864A0','NMPC':'#CE691C','A_NMPC':'#218459'}
STYLES={'ARC':'-','NMPC':(0,(6,3)),'A_NMPC':(0,(6,2,1.5,2))}
NAMES={'ARC':'ARC','NMPC':'NMPC','A_NMPC':'A-NMPC'}

def style(size=10):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':size,'axes.labelsize':size,
       'axes.titlesize':size,'xtick.labelsize':size-1,'ytick.labelsize':size-1,
       'legend.fontsize':size-1,'axes.spines.top':False,'axes.spines.right':False,
       'pdf.fonttype':42,'ps.fonttype':42,'savefig.bbox':'tight'})

def save(fig,name):
    fig.savefig(D/(name+'.pdf'));fig.savefig(D/(name+'.png'),dpi=160);plt.close(fig)

def trajectory(scenarios,name):
    columns=len(scenarios);scale=columns
    style(10*scale)
    fig,axes=plt.subplots(5,columns,figsize=(5.4*scale,6.2*scale),squeeze=False,layout='constrained')
    labels=[]
    for j,scenario in enumerate(scenarios):
        maximum=0;p_peak=1.6;p_min=.1;pressure_detail_end=3300/3600
        for controller in ['ARC','NMPC','A_NMPC']:
            row,data=load_run(run_folder(scenario,controller,200))
            t=data['t']/3600;x=data['x'];maximum=max(maximum,t[-1]);lw=1.5*scale
            p_peak=max(p_peak,float(np.max(x[:,13])/1e6));p_min=min(p_min,float(np.min(x[:,13])/1e6)-.02)
            for k in [0,1]:
                line=axes[k,j].plot(t,x[:,11],color=COLORS[controller],ls=STYLES[controller],lw=lw,label=NAMES[controller])[0]
            axes[2,j].plot(t,x[:,13]/1e6,color=COLORS[controller],ls=STYLES[controller],lw=lw)
            mask=(data['t']>=3300)&(x[:,15]<3250-1e-3)
            if np.any(mask):
                pressure_detail_end=max(pressure_detail_end,float(t[mask][-1]))
                axes[3,j].plot(t[mask],x[mask,13]/1000-1500,color=COLORS[controller],ls=STYLES[controller],lw=lw)
            # 30 s block means retain the same visual resolution for all controllers.
            ut=t[:-1];uf=data['u'][:,1]/4e-7;n=len(ut)
            at=np.array([ut[k] for k in range(0,n,30)])
            au=np.array([np.mean(uf[k:k+30]) for k in range(0,n,30)])
            axes[4,j].step(at,au,where='post',color=COLORS[controller],ls=STYLES[controller],lw=lw)
            if controller=='ARC':
                axes[2,j].plot(t[:-1],data['Psp']/1e6,color='#873B99',ls=(0,(1.2,2)),lw=1.7*scale,label='ARC pressure setpoint')
                axes[3,j].plot(t[:-1][mask[:-1]],data['Psp'][mask[:-1]]/1000-1500,color='#873B99',ls=(0,(1.2,2)),lw=1.7*scale)
            if not row['completed']:
                axes[0,j].plot(t[-1],x[-1,11],'x',color=COLORS[controller],ms=6*scale,mew=1.3*scale)
        for k in [0,1]:
            axes[k,j].axhline(351.7,c='black',ls=(0,(4,3)),lw=1.0*scale,zorder=0)
            axes[k,j].axhline(351.,c='#888888',ls=':',lw=1.0*scale,zorder=0)
        axes[1,j].set(ylim=(349.8,352.4),xlim=(3300/3600,maximum+.03))
        axes[1,j].axhspan(350.3,351.7,color='#777777',alpha=.09,zorder=-1)
        axes[2,j].axhline(1.6,c='black',ls=(0,(4,3)),lw=1.0*scale,zorder=0)
        axes[2,j].set(ylim=(p_min,max(1.65,p_peak+.03)))
        axes[3,j].set_xlim(3300/3600,pressure_detail_end+.03)
        axes[3,j].axhline(0,c='#888888',ls=':',lw=.8*scale,zorder=0)
        axes[4,j].set(ylim=(-.04,1.1),xlabel='Time (h)')
        axes[0,j].set_xlim(0,maximum+.03)
        axes[2,j].set_xlim(0,maximum+.03);axes[4,j].set_xlim(0,maximum+.03)
        title={'N':'N: nominal','PM_plus':'PM+: adverse mismatch','F':'F: mismatch + gel effect'}[scenario]
        for k,lab in enumerate(['Reactor T (K)','T detail (K)','Pressure (MPa)',r'$\Delta P$ (kPa)',r'$u_B/u_B^{\max}$']):
            axes[k,j].set_ylabel(lab)
            axes[k,j].set_title('('+string.ascii_lowercase[k*columns+j]+') '+(title if k==0 else ['','Operating limit detail','Pressure and setpoint','Pressure detail','Applied initiator'][k]),loc='left')
            axes[k,j].grid(alpha=.22)
        if columns==1:axes[2,j].legend(loc='lower right')
    handles,labels=axes[0,0].get_legend_handles_labels()
    if columns==2:
        from matplotlib.lines import Line2D
        handles.append(Line2D([0],[0],color='#873B99',ls=(0,(1.2,2)),lw=3));labels.append('ARC pressure setpoint')
    fig.legend(handles,labels,loc='outside lower center',ncol=2 if columns==2 else 3,frameon=False)
    save(fig,name)

def signals():
    style(10)
    fig,ax=plt.subplots(3,2,figsize=(5.4,6.5),layout='constrained')
    for j,scenario in enumerate(['N','PM_plus','F']):
        row,data=load_run(OUT/f'{scenario}_ARC_seed000');t=data['t'][:-1]/3600
        ax[j,0].plot(t,data['u'][:,2],color=COLORS['ARC'],lw=1.3)
        ax[j,1].plot(t,data['virtual'],color='#873B99',lw=1.3)
        for k in range(2):
            ax[j,k].axhline(-1,c='black',ls='--',lw=1)
            ax[j,k].set_title(f'({string.ascii_lowercase[j*2+k]}) '+scenario.replace('PM_plus','PM+'),loc='left')
            ax[j,k].grid(alpha=.2)
        ax[j,0].set_ylabel('Applied utility input');ax[j,0].set_ylim(-1.1,1.1)
        ax[j,1].set_ylabel('Virtual split-range demand')
        ax[j,1].axhline(-.95,c='#777777',ls=':',lw=1)
    for k in range(2):ax[-1,k].set_xlabel('Time (h)')
    save(fig,'industrial_utility_signals')

def diagrams():
    """Compile local corrections to the submitted process and recipe drawings."""
    import subprocess
    for name in ['control_structure_revision','mode_transition_revision']:
        subprocess.run(['latexmk','-norc','-pdf','-interaction=nonstopmode','-halt-on-error','-outdir='+str(output_dir()),name+'.tex'],cwd=ROOT/'docs/figure_sources',check=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('part',choices=['diagrams','trajectories','signals']);a=p.parse_args()
    if a.part=='diagrams':diagrams()
    elif a.part=='signals':signals()
    else:
        trajectory(['PM_plus'],'industrial_pmplus_revision')
        trajectory(['N','F'],'industrial_nominal_fault_revision');signals()
