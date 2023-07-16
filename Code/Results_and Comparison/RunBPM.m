function RunBPM(dir_name,in_data_file,n_subd, name_ext)
% function RunBPM(dir_name,in_off_file,out_off_file, show_figs_flag, save_figs_flag,save_off_flag,comparison_flag, cauchy_flag)
tdata = load('planarity_colormap.mat');

[~,f_in,~] = fileparts(in_data_file.name);
[Vz,Vw,F] = ExtractMeshesFromUV(f_in);
mesh_name = erase(f_in,name_ext);

% mesh_name = erase(f_in,"_in");
Mz = []; Mz.V = Vz; Mz.F = F;
Mw = []; Mw.V = Vw; Mw.F = F;

% new_dir = [dir_name,'data_out/',mesh_name,'/'];
new_dir = [dir_name,mesh_name,'/'];

if ~exist(new_dir, 'dir')
    mkdir(new_dir);
end
% copy the texture file to the folder
I = imread('texture.png');
imwrite(I,[new_dir,'texture.png'])
%% set parameters

limits_qc = [1 1.1]; % for KVF limits_qc=[1 2]


%% generate colormap figure

cbarfigname = [new_dir 'colorbar'];
fh = figure('visible','off');
colormap(fh,tdata.cm1);
caxis(limits_qc);
set(fh,'color','w'); axis off
c = colorbar; c.FontSize = 15; c.TickLabelInterpreter='latex';
export_fig(cbarfigname ,'-m3');


%% Mz 
filename = [new_dir,mesh_name,'_Mz','.png'];
SaveMeshFigs(Mz, filename, true, false);

%% Mw

filename = [new_dir,mesh_name,'_Mw','.png'];
SaveMeshFigs(Mw, filename, true, false);


%% QC error between Mz and Mw
[qc_err_Mz2Mw, ~, max_qc, avg_qc] = ComputeQuasiConformalError(Mz.V,Mw.V,Mz.F);
    disp([mesh_name,'   Orig:   max=', num2str(max_qc),'   avg_qc=',num2str(avg_qc)])
qc_err_Mz2Mw = real(qc_err_Mz2Mw);
filename = [new_dir,mesh_name,'_Mz2Mw_Mw_QCerr'];
SaveQCErrorFigs(Mw, qc_err_Mz2Mw,limits_qc,tdata.cm1 ,filename, true, false);

filename = [new_dir,mesh_name,'_Mz2Mw_Mz_QCerr'];
SaveQCErrorFigs(Mz, qc_err_Mz2Mw,limits_qc,tdata.cm1 ,filename, true, false);

%% Compute piecewise continuous conformal mapping - BPM
disp('running BPMContinuousParam')
[Mw_fine_BPM, Mz_fine_BPM] = BPMContinuousParam(Vz, Vw, F ,n_subd);


%% Compute piecewise continuous conformal mapping - Projective
disp('running ProjectiveContinuousParam')
[Mw_fine_proj, Mz_fine_proj] = ProjectiveContinuousParam(Vz, Vw, F ,n_subd);

%% Compute piecewise continuous conformal mapping - Linear
disp('running BPMContinuousParam')
[Mw_fine_linear, Mz_fine_linear] = LinearContinuousParam(Vz, Vw, F ,n_subd);
%% texture mapping 
texture_filename = 'texture.png';
ext = '.obj';

%% BPM - texture mapping 

% on the source

name = 'Mz_BPM' ;
M.vertices = Mz_fine_BPM.V; M.triangles = Mz_fine_BPM.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% UV_w = NormalizeUV(Mw_fine_BPM.V(:,1:2),Mw.V(:,1:2));
UV_w = NormalizeUV(Mw_fine_BPM.V(:,1:2),Mw_fine_BPM.V(:,1:2));
generate_texture_mesh(M, UV_w, texture_filename, output_M_filename)

UV_z = NormalizeUV(Mz_fine_BPM.V(:,1:2),Mw_fine_BPM.V(:,1:2));

% name = 'Mz_texture';
% output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% % UV_z = NormalizeUV(Mz_fine_BPM.V(:,1:2),Mz.V(:,1:2));
% % UV_z = NormalizeUV(Mz_fine_BPM.V(:,1:2));
% generate_texture_mesh(M, UV_z, texture_filename,output_M_filename)

% on the target
name = 'Mw_BPM';
M.vertices = Mw_fine_BPM.V; M.triangles = Mw_fine_BPM.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
generate_texture_mesh(M, UV_z, texture_filename,output_M_filename)

% % save the texture mesh 
% name = 'Mw_texture';
% output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% generate_texture_mesh(M, UV_w, texture_filename,output_M_filename)

%% Show texture mapping -Projective

% Proj - show the mapping on the source

name = 'Mz_Proj' ;
M.vertices = Mz_fine_proj.V; M.triangles = Mz_fine_proj.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% UV = NormalizeUV(Mw_fine_proj.V(:,1:2));
% UV = NormalizeUV(Mw_fine_proj.V(:,1:2));
UV = NormalizeUV(Mw_fine_proj.V(:,1:2),Mw_fine_BPM.V(:,1:2));

generate_texture_mesh(M, UV, texture_filename, output_M_filename)

% show the mapping on the target
name = 'Mw_Proj';
M.vertices = Mw_fine_proj.V; M.triangles = Mw_fine_proj.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% UV = NormalizeUV(Mz_fine_proj.V(:,1:2));
UV = NormalizeUV(Mz_fine_proj.V(:,1:2),Mw_fine_BPM.V(:,1:2));

generate_texture_mesh(M, UV, texture_filename,output_M_filename)

%% Show texture mapping - Linear

% linear - show the mapping on the source
name = 'Mz_Linear' ;
M.vertices = Mz_fine_linear.V; M.triangles = Mz_fine_linear.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% UV = NormalizeUV(Mw_fine_linear.V(:,1:2));
UV = NormalizeUV(Mw_fine_linear.V(:,1:2),Mw_fine_BPM.V(:,1:2));

generate_texture_mesh(M, UV, texture_filename, output_M_filename)

% show the mapping on the target
name = 'Mw_Linear';
M.vertices = Mw_fine_linear.V; M.triangles = Mw_fine_linear.F;
output_M_filename = [new_dir,mesh_name,'_',name,'_out',ext];
% UV = NormalizeUV(Mz_fine_linear.V(:,1:2));
UV = NormalizeUV(Mz_fine_linear.V(:,1:2),Mw_fine_BPM.V(:,1:2));

generate_texture_mesh(M, UV, texture_filename,output_M_filename)






%% quasi-conformal distortion - BPM

 [qc_error, ~, max_qc, avg_qc] = ComputeQuasiConformalError(Mz_fine_BPM.V,Mw_fine_BPM.V,Mz_fine_BPM.F);
 qc_error = real(qc_error);
 disp([mesh_name,'   BPM:   max=', num2str(max_qc),'   avg_qc=',num2str(avg_qc)])
 filename = [new_dir,mesh_name,'_Mz2Mw_Mw_QCerr_BPM'];
 SaveQCErrorFigs(Mw_fine_BPM, qc_error,limits_qc,tdata.cm1 ,filename, true, false);

 
 filename = [new_dir,mesh_name,'_Mz2Mw_Mz_QCerr_BPM'];
 SaveQCErrorFigs(Mz_fine_BPM, qc_error,limits_qc,tdata.cm1 ,filename, true, false);


%% quasi-conformal distortion - Projective

[qc_error_proj, ~, max_qc, avg_qc] = ComputeQuasiConformalError(Mz_fine_proj.V,Mw_fine_proj.V,Mz_fine_proj.F);
 qc_error_proj = real(qc_error_proj);

filename = [new_dir,mesh_name,'_Mz2Mw_Mw_QCerror_proj'];
SaveQCErrorFigs(Mw_fine_proj, qc_error_proj,limits_qc,tdata.cm1 ,filename, true, false);

filename = [new_dir,mesh_name,'_Mz2Mw_Mz_QCerror_proj'];
SaveQCErrorFigs(Mz_fine_proj, qc_error_proj,limits_qc,tdata.cm1 ,filename, true, false);


%% quasi-conformal distortion - Linear

[qc_error_linear, ~, max_qc, avg_qc] = ComputeQuasiConformalError(Mz_fine_linear.V,Mw_fine_linear.V,Mz_fine_linear.F);
% SaveQCErrorFigs(Mw_fine_proj, qc_error_proj,limits_qc,tdata.cm1 ,filename, save_figs_flag, show_figs_flag);
qc_error_linear = real(qc_error_linear);
filename = [new_dir,mesh_name,'_Mz2Mw_Mw_QCerror_linear'];
SaveQCErrorFigs(Mw_fine_linear, qc_error_linear,limits_qc,tdata.cm1 ,filename, true, false);

filename = [new_dir,mesh_name,'_Mz2Mw_Mz_QCerror_linear'];
SaveQCErrorFigs(Mz_fine_linear, qc_error_linear,limits_qc,tdata.cm1 ,filename, true, false);




end