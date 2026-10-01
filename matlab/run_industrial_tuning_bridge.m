function run_industrial_tuning_bridge(root)
% Independent prospective commissioning study; frozen controller is reused.
addpath(fullfile(root,'.runtime','casadi','matlab'));
addpath(fullfile(root,'matlab'));
protocol=jsondecode(fileread(fullfile(root,'state','industrial_tuning_bridge.json')));
seed=protocol.fresh_batches.seed;
ndata=load(fullfile(root,'results','jpc_revision_20260907',sprintf('noise_%03d.mat',seed)),'paired_noise');
cases=protocol.fresh_batches.candidates_K_Ti;
scenarios=protocol.fresh_batches.scenarios;
for j=1:numel(scenarios)
    scenario=scenarios{j};
    for k=1:size(cases,1)
        controller=sprintf('ARC_K%d_T%d',cases(k,1),cases(k,2));
        folder=fullfile(root,'results','jpc_revision_20260907','industrial_tuning_bridge',sprintf('%s_%s_seed%03d',scenario,controller,seed));
        if exist(fullfile(folder,'trajectory.mat'),'file'), error('Refusing to overwrite %s',folder); end
        fault=strcmp(scenario,'F');
        [sys,par]=ptfe(struct('tf',1,'isFault',fault));
        dk=par.d0;
        if strcmp(scenario,'PM_plus') || fault, dk(1:2)=dk(1:2).*[1.2;.8]; end
        par.seed=seed; par.max_steps=protocol.fresh_batches.max_steps;
        par.noise_scale=1; par.arc_variant='full_arc';
        par.initial_actuator=[.694;0;0]; par.paired_noise=ndata.paired_noise;
        par.plot_fig=false; par.K_VPC=cases(k,1); par.Ti_VPC=cases(k,2);
        clock_start=tic;
        [~,~,SimData]=arc_pid_revision(sys,par,dk);
        elapsed=toc(clock_start);
        if ~exist(folder,'dir'), mkdir(folder); end
        save(fullfile(folder,'trajectory.mat'),'SimData','dk','elapsed','seed','scenario','controller','protocol','-v7');
        fprintf('TUNING %s %s: %.4f h, peak %.4f K, %.2f s wall\n',controller,scenario,size(SimData.X,2)/3600,max(SimData.X(12,:)),elapsed);
    end
end
end
