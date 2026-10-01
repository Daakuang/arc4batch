"""Temperature-quality constrained NMPC; pressure tracking has secondary priority."""
import time
import numpy as np
import casadi as ca
from adaptive_model import SCALE_U,LOW_U,INIT_U,arr

class QualityController:
    def __init__(self,model,horizon=60):
        self.model=model;self.N=horizon;self.last=INIT_U/SCALE_U
        roots=np.r_[0.,ca.collocation_points(3,'radau')]
        C=np.zeros((4,4));D=np.zeros(4);B=np.zeros(4)
        for j in range(4):
            poly=np.poly1d([1.])
            for r in range(4):
                if r!=j:poly*=np.poly1d([1.,-roots[r]])/(roots[j]-roots[r])
            D[j]=poly(1.);C[j,:]=np.polyder(poly)(roots)
            ip=np.polyint(poly);B[j]=ip(1.)-ip(0.)
        p=ca.MX.sym('parameters',19)
        variables=[];lower=[];upper=[];guess=[];g=[];gl=[];gu=[]
        def variable(name,n,lo,hi,start):
            q=ca.MX.sym(name,n);variables.append(q)
            lower.extend(np.broadcast_to(lo,(n,)));upper.extend(np.broadcast_to(hi,(n,)))
            guess.extend(np.broadcast_to(start,(n,)));return q
        def constrain(q,lo=0.,hi=0.):
            n=int(q.numel());g.append(q);gl.extend(np.broadcast_to(lo,(n,)));gu.extend(np.broadcast_to(hi,(n,)))
        eps=variable('operating_slacks',2,[0.,0.],[21.45,.4],[.001,.00001])
        x=variable('initial_state',12,-np.inf,np.inf,model.initial)
        constrain(x-p[:12])
        self.u_indices=[];self.state_indices=[2]
        objective=p[17]*(1e5*(eps[0]+eps[0]**2)+1e7*(eps[1]+eps[1]**2))
        previous=p[12:15]
        def operating(q,future_s):
            constrain(ca.vertcat(q[5]*model.scale[5]-351.7-eps[0],
                       ca.if_else(p[18]+future_s>=3300.,350.3,320.15)-q[5]*model.scale[5]-eps[0],
                       (q[7]*model.scale[7]-1.6e6)/1e6-eps[1]),-np.inf,0.)
        for k in range(horizon):
            self.u_indices.append(len(lower))
            u=variable(f'u{k}',3,LOW_U/SCALE_U,1.,INIT_U/SCALE_U)
            stages=[]
            for j in range(3):
                self.state_indices.append(len(lower))
                z=variable(f'x{k}_{j}',12,model.lo,model.hi,model.initial)
                stages.append(z);operating(z,30*(k+roots[j+1]))
            end=D[0]*x
            for j in range(1,4):
                derivative=C[0,j]*x
                for r in range(1,4):derivative+=C[r,j]*stages[r-1]
                dx,cost=model.f(stages[j-1],u,p[15:17])
                cost+=ca.if_else(p[18]+30*(k+roots[j])>=3300.,model.hold_penalty(stages[j-1]),0.)
                constrain(30*dx-derivative)
                objective+=30*B[j]*cost;end+=D[j]*stages[j-1]
            objective+=ca.dot(ca.DM([8,.1,4]),((previous-u)*SCALE_U)**2)/1000.
            previous=u
            self.state_indices.append(len(lower))
            x=variable(f'end{k}',12,model.lo,model.hi,model.initial)
            constrain(end-x);operating(x,30*(k+1))
        w=ca.vertcat(*variables);g=ca.vertcat(*g)
        options=dict(print_time=False,error_on_fail=False,expand=True,
            ipopt=dict(print_level=0,sb='yes',max_iter=500,tol=1e-6,acceptable_tol=1e-5,
                       acceptable_iter=5,bound_relax_factor=1e-9))
        self.solver=ca.nlpsol('parameter_nmpc','ipopt',dict(x=w,p=p,f=objective,g=g),options)
        self.constraint=ca.Function('parameter_constraints',[w,p],[g])
        self.objective=ca.Function('quality_objective',[w,p],[objective])
        self.lb=np.array(lower);self.ub=np.array(upper);self.w=np.array(guess)
        self.gl=np.array(gl);self.gu=np.array(gu)

    def control(self,estimate,t_s=0):
        m=self.model
        for channel,state,cap,tol in [(0,9,3250.,1e-3),(1,10,.002,1e-10)]:
            if estimate[state]*m.scale[state]>=cap-tol:
                for index in self.u_indices:self.ub[index+channel]=0.;self.w[index+channel]=0.
        already_above=bool(estimate[5]*m.scale[5]>351.7 or estimate[7]*m.scale[7]>1.6e6 or (t_s>=3300 and estimate[5]*m.scale[5]<350.3))
        p=np.r_[estimate[:12],self.last,estimate[12:],float(already_above),float(t_s)]
        self.w[2:14]=estimate[:12]
        attempts=[];start=time.perf_counter();accepted=False;solution=None
        for attempt in range(2):
            recovery=already_above or attempt>0
            p[17]=float(recovery)
            self.ub[:2]=[21.45,.4] if recovery else [0.,0.]
            if not recovery:self.w[:2]=0.
            sol=self.solver(x0=self.w,p=p,lbx=self.lb,ubx=self.ub,lbg=self.gl,ubg=self.gu)
            w=arr(sol['x']);stats=self.solver.stats();gg=arr(self.constraint(w,p))
            violation=max(0.,float(np.max(self.gl-gg)),float(np.max(gg-self.gu)),
                          float(np.max(self.lb-w)),float(np.max(w-self.ub)))
            accepted=bool(stats['success'] and np.isfinite(w).all() and violation<1e-6)
            info=dict(status=stats['return_status'],mode='recovery' if recovery else 'hard',iterations=int(stats['iter_count']),
                      residual=violation,accepted=accepted,objective=float(sol['f']),
                      temperature_slack_K=float(max(0,w[0])),pressure_slack_MPa=float(max(0,w[1])))
            attempts.append(info)
            if accepted:solution=w;break
            # Reinitialize a failed NLP using the current state at every node.
            # This retry is independent of the failed candidate's feed request.
            self.w[0]=max(.001,estimate[5]*m.scale[5]-351.7+.01,350.3-estimate[5]*m.scale[5]+.01 if t_s>=3300 else 0.)
            self.w[1]=max(.00001,(estimate[7]*m.scale[7]-1.6e6)/1e6+.001)
            for i in self.state_indices:self.w[i:i+12]=estimate[:12]
            for i in self.u_indices:self.w[i:i+3]=np.minimum(self.ub[i:i+3],self.last)
        if accepted:
            self.last=np.clip(solution[self.u_indices[0]:self.u_indices[0]+3],LOW_U/SCALE_U,1.)
            block=51
            self.w=np.r_[solution[:2],solution[2+block:14+block],solution[14+block:],solution[-block:]]
        else:
            self.last=np.array([0.,0.,-1.])
        info=dict(attempts[-1],seconds=time.perf_counter()-start,attempts=attempts,
                  fallback='none' if accepted else 'feed_stop_max_cooling')
        return self.last*SCALE_U,info
