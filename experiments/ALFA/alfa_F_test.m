% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
% This is new!
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;

addpath('surfstat/')
addpath('utils/')Image
set(groot,'defaultFigureVisible','off')

%% folders where meshes are located
meshes_l = "";
meshes_r = "";

%Csv with the information of the meshes
out_csv = "";

% Path of the average, in .obj format
avg_l_path = '';
avg_r_path = '';
% Path of the average, in .vtk format
avg_l_vtk = '';
avg_r_vtk = '';

% Base exp to store the results
exp_dir_base = '';
mkdir(exp_dir_base);

data = readtable(out_csv);
% Number of subjects
N = size(data,1);

%% Load the data
Data_loading;

%% Models
% Left hippocampus
% Create and compute model
% Both for left and right hippocampus

t_l = double(mesh_l.coord);
t_lavg = double(avg_l.coord)';
t_lavg=kron(t_lavg,ones(N,1));
t_lavg=reshape(t_lavg,N,2511,3);
t_l = t_l - t_lavg;

t_r = double(mesh_r.coord);
t_ravg = double(avg_r.coord)';
t_ravg=kron(t_ravg,ones(N,1));
t_ravg=reshape(t_ravg,N,2636,3);
t_r = t_r - t_ravg;

% Baseline model
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe;
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);


% Interaction model
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe + Age*Apoe + Agesq*Apoe;
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe + Age*Apoe + Agesq*Apoe;

slm_l_int = SurfStatLinMod(t_l, model_l, avg_l);
slm_r_int = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int + Age*Apoe_int + Agesq*Apoe_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int + Age*Apoe_int + Agesq*Apoe_int;

slm_l_add_int = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add_int = SurfStatLinMod(t_r, model_r_add, avg_r);


%% F TESTING
% For each test, save to disk the original map and the FDR corrected map.
% I will end up writing a function for this... maybe lmao

models_base = {slm_r, slm_r_add_int, };

model_extended = {slm_r_int, slm_r_add, };

avg = {avg_r, avg_r, };

title = {'F statistic, right hipp, FDR', 'F statistic, additive right hipp'};

out_dir = {'fstat_r', 'fstat_r_fdr', 'fstat_add_r' };

for i = 1:length(contrast_list)
 % also create the .vtk
F_stats = SurfStatF( models_base{i}, model_extended{i} );
SurfStatView( F_stats.t(1,:), avg{i}, title{i});
saveas(gcf,strcat(exp_dir_base, out_dir{i}, '.png'))
close

SurfStatView( SurfStatQ( slm_add_r, [] ), avg_r, 'F statistic, additive right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_fdr.png'))
close

csv_out = strcat(exp_dir, save_file{i},'_values.csv');
to_csv_mat = [F_stats.t(1,:)', F_stats.t(1,:)', F_stats.t'];
csvwrite(csv_out, to_csv_mat);

% Call the function from python to create the volume
vtk_out = strcat(exp_dir, save_file{i},'_Tvalues.vtk');
system("python utils/create_vtk_from_matlab.py " + avg_vtk + " " +...
       csv_out + " " + vtk_out);
end

slm_add_r = SurfStatF( slm_r_add_int, slm_r_add );
SurfStatView( slm_add_r.t(1,:), avg_r, 'F statistic, additive right hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r.png'))
close

SurfStatView( SurfStatQ( slm_add_r, [] ), avg_r, 'F statistic, additive right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_fdr.png'))
close

slm_base_l = SurfStatF( slm_l_int, slm_l );
SurfStatView( slm_base_l.t(1,:), avg_l, 'F statistic, left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l.png'))
close

SurfStatView( SurfStatQ( slm_base_l, [] ), avg_l, 'F statistic, left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_fdr.png'))
close

slm_add_l = SurfStatF( slm_l_add_int, slm_l_add );
SurfStatView( slm_add_l.t(1,:), avg_l, 'F statistic, additive left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l.png'))
close

SurfStatView( SurfStatQ( slm_add_l, [] ), avg_l, 'F statistic, additive left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_fdr.png'))
close
