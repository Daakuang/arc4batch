function export_models(root)
import casadi.*
addpath(fullfile(root,'.runtime','casadi','matlab'));
addpath(fullfile(root,'matlab'));
dest=fullfile(root,'results','jpc_revision_20260907','models');
if ~exist(dest,'dir'), mkdir(dest); end
for fault=0:1
    [s,p]=ptfe(struct('tf',1,'isFault',fault));
    s.F.save(fullfile(dest,sprintf('plant_flow_%d.casadi',fault)));
    s.m.save(fullfile(dest,sprintf('plant_m_%d.casadi',fault)));
    s.f.save(fullfile(dest,sprintf('plant_rhs_%d.casadi',fault)));
    data=struct('x0',p.x0,'u0',p.u0,'d0',p.d0,'lbu',p.lbu,'ubu',p.ubu);
    fid=fopen(fullfile(dest,sprintf('plant_%d.json',fault)),'w');
    fprintf(fid,'%s',jsonencode(data)); fclose(fid);
end
[s,p]=ptfe_reduced(struct('tf',1,'isFault',0,'SMS',1));
s.f.save(fullfile(dest,'reduced_f.casadi'));
data=struct('x0',p.x0,'u0',p.u0,'d0',p.d0,'lbu',p.lbu,'ubu',p.ubu, ...
    'lbx',p.lbx,'ubx',p.ubx,'x_norm',p.x_norm,'u_norm',p.u_norm);
fid=fopen(fullfile(dest,'reduced.json'),'w'); fprintf(fid,'%s',jsonencode(data)); fclose(fid);
fprintf('Exported exact CasADi models, version %s\n',CasadiMeta.version());
end
