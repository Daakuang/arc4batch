function run_arc_revision(root,scenario,seed,max_steps,variant)
addpath(fullfile(root,'.runtime','casadi','matlab'));
addpath(fullfile(root,'matlab'));
if nargin<4, max_steps=30000; end
if nargin<5, variant='full_arc'; end
fault=strcmp(scenario,'F');
[sys,par]=ptfe(struct('tf',1,'isFault',fault));
dk=par.d0;
if strcmp(scenario,'PM_plus') || fault, dk(1:2)=dk(1:2).*[1.2;.8]; end
if strcmp(scenario,'PM_minus'), dk(1:2)=dk(1:2).*[.8;1.2]; end
par.seed=seed; par.max_steps=max_steps; par.noise_scale=1;
par.arc_variant=variant;
par.initial_actuator=[.694;0;0];
noise_file=fullfile(root,'results','jpc_revision_20260907',sprintf('noise_%03d.mat',seed));
ndata=load(noise_file,'paired_noise'); par.paired_noise=ndata.paired_noise;
par.plot_fig=false;
clock_start=tic;
[~,~,SimData]=arc_pid_revision(sys,par,dk);
elapsed=toc(clock_start);
controller='ARC';
if ~strcmp(variant,'full_arc'), controller=['ARC_' variant]; end
folder=fullfile(root,'results','jpc_revision_20260907',sprintf('%s_%s_seed%03d',scenario,controller,seed));
if ~exist(folder,'dir'), mkdir(folder); end
save(fullfile(folder,'trajectory.mat'),'SimData','dk','elapsed','seed','scenario','variant','-v7');
fprintf('%s %s seed %d: %.4f h, peak %.4f K, CPU %.2f s\n',controller,scenario,seed,size(SimData.X,2)/3600,max(SimData.X(12,:)),elapsed);
end
