function run_arc_additional(root)
addpath(fullfile(root,'.runtime','casadi','matlab'));
addpath(fullfile(root,'matlab'));
base=fullfile(root,'results','study');
for j=1:5
    [sys,par]=ptfe(struct('tf',1,'isFault',false));
    par.max_steps=30000; par.initial_actuator=[.694;0;0]; par.plot_fig=false;
    par.arc_variant='full_arc'; dk=par.d0;
    if j<=3
        scenario='PM_plus';seed=94;kind='tuning';names={'conservative','baseline','aggressive'};
        gains=[6,8,10];targets=[-.90,-.95,-.98];
        par.K_VPC=gains(j);par.VPC_Setpoint=targets(j);par.noise_scale=1;
        controller=['ARC_' names{j}];
    else
        names={'N','PM_plus'};scenario=names{j-3};seed=50;kind='stress';controller='ARC';
        par.noise_scale=3;
    end
    if strcmp(scenario,'PM_plus'),dk(1:2)=dk(1:2).*[1.2;.8];end
    par.seed=seed;nd=load(fullfile(base,sprintf('noise_%03d.mat',seed)),'paired_noise');
    par.paired_noise=nd.paired_noise;clock=tic;
    if j<=3,[~,~,SimData]=arc_controller(sys,par,dk);
    else,[~,~,SimData]=arc_controller_stress(sys,par,dk);end
    elapsed=toc(clock);
    folder=fullfile(base,kind,sprintf('%s_%s_seed%03d',scenario,controller,seed));
    if ~exist(folder,'dir'),mkdir(folder);end
    variant=controller;
    save(fullfile(folder,'trajectory.mat'),'SimData','dk','elapsed','seed','scenario','variant','-v7');
    fprintf('%s %s seed %d: %.4f h, peak %.4f K, CPU %.2f s\n',controller,scenario,seed,size(SimData.X,2)/3600,max(SimData.X(12,:)),elapsed);
end
end
