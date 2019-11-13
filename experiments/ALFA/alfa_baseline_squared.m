% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;
clear all
close all
addpath('surfstat/')
addpath('utils/')
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

model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe;
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);

% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and volume.
corr = [1 4 5 6];
corr_int = [1 4 5 6];

%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)

%% Left

contrast_list = {-age, age, -agesq, agesq, volume_l, -volume_l,...
                 Gender.M-Gender.F, Gender.F-Gender.M,...
                 Apoe.HE-Apoe.NC, Apoe.NC-Apoe.HE,...
                 Apoe.HO-Apoe.NC, Apoe.NC-Apoe.HO,...
                 Apoe.HE-Apoe.HO, Apoe.HO-Apoe.HE,...
                 (0.5*Apoe.HE+0.5*Apoe.NC) - Apoe.HO,...
                 Apoe.HO - (0.5*Apoe.HE+0.5*Apoe.NC),...
                 (0.5*Apoe.HE+0.5*Apoe.HO) - Apoe.NC,...
                 Apoe.NC - (0.5*Apoe.HE+0.5*Apoe.HO)};
    
% List of files to save
save_file = {'-age', '+age', '-agesq', 'agesq', '+volume', '-volume',...
             '+male_-female', '+female_-male',... % '-yed', 'yed',...
             '+HE_-NC', '+NC_-HE', '+HO_-NC', '+NC_-HO', '+HO_-HE',...  
             '+HE_-HO', 'HE+NC-HO', '-HE-NC+HO', 'HE+HO-NC', '-HE-HO+NC'};

% ADDITIVE MODEL
contrast_list_add = {apoe_int, -apoe_int};
save_file_add = {'NC_HE_HO', '-NC_HE_HO'};

% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo/');
mkdir(exp_dir);

% Save model of age + age^2 to disk
% hard coded, will move to a function later

csv_out = strcat(exp_dir, 'agesq_values.csv');
to_csv_mat = [squeeze(slm_l.coef(2,:,:)), squeeze(slm_l.coef(3,:,:))];
csvwrite(csv_out, to_csv_mat);

% Call the function from python to create the volume
vtk_out = strcat(exp_dir,'agesq_effects.vtk');
system("python utils/create_squared_vtk.py " + avg_l_vtk + " " +...
       csv_out + " " + vtk_out);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, apoe, corr)

% ADDITIVE MODEL
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, age, apoe, corr_int)

 


%% Right


% List of contrasts
contrast_list = {-age, age, -agesq, agesq, volume_r, -volume_r,...
                 Gender.M-Gender.F, Gender.F-Gender.M,...
                 Apoe.HE-Apoe.NC, Apoe.NC-Apoe.HE,...
                 Apoe.HO-Apoe.NC, Apoe.NC-Apoe.HO,...
                 Apoe.HE-Apoe.HO, Apoe.HO-Apoe.HE,...
                 (0.5*Apoe.HE+0.5*Apoe.NC) - Apoe.HO,...
                 Apoe.HO - (0.5*Apoe.HE+0.5*Apoe.NC),...
                 (0.5*Apoe.HE+0.5*Apoe.HO) - Apoe.NC,...
                 Apoe.NC - (0.5*Apoe.HE+0.5*Apoe.HO)};


% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

csv_out = strcat(exp_dir, 'agesq_values.csv');
to_csv_mat = [squeeze(slm_r.coef(2,:,:)), squeeze(slm_r.coef(3,:,:))];
csvwrite(csv_out, to_csv_mat);

% Call the function from python to create the volume
vtk_out = strcat(exp_dir,'agesq_effects.vtk');
system("python utils/create_squared_vtk.py " + avg_r_vtk + " " +...
       csv_out + " " + vtk_out);

% Call save function
save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, apoe, corr)


% additive
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, apoe, corr_int)


%% INTERACTIONS TIME
