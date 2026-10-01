function run_arc_suite(root,seeds,ablation)
if nargin<3, ablation=false; end
addpath(fullfile(root,'matlab'));
scenarios={'N','PM_minus','PM_plus','F'};
variants={'full_arc'};
if ablation
    scenarios={'N','PM_plus','F'};
    variants={'fixed_psp','recipe_fb','recipe_only'};
end
for seed=seeds
    for j=1:numel(scenarios)
        for k=1:numel(variants)
            simulate_arc(root,scenarios{j},seed,30000,variants{k});
        end
    end
end
end
