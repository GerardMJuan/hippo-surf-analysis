% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;
clear all
close all
addpath('../../surfstat/')
addpath('../../utils/')
set(groot,'defaultFigureVisible','off')

%% folders where meshes are located
meshes_l = "";
meshes_r = "";

%Csv with the information of the meshes
csv = "";

% Path of the average, in .obj format
avg_l_path = '';
avg_r_path = '';
% Path of the average, in .vtk format
avg_l_vtk = '';
avg_r_vtk = '';

% Base exp to store the results
exp_dir_base = '';
mkdir(exp_dir_base);

data = readtable(csv);
% Number of subjects
N = size(data,1);

%% Load covariates
Data_loading_all;

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

model_l = 1 + Age + Gender + Volume_l + Site + Apoe + DX;
model_r = 1 + Age + Gender + Volume_r + Site + Apoe + DX;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Gender + Volume_r + Site + Apoe_int + DX_int;
model_l_add = 1 + Age + Gender + Volume_l + Site + Apoe_int + DX_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);

% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and 
% the three diagnosis.

% for the additive part, is different, as we only have a covariate.

% CONFIRMAR
corr = [];
corr_int = [];

%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)


% NO INTERACTIONS
% List of contrasts
contrast_list = {-age, age, site, -site...
                 Gender.M-Gender.F, Gender.F-Gender.M,...
                 -volume_l, volume_l,...
                 Apoe.HE-Apoe.NC, Apoe.NC-Apoe.HE,...
                 Apoe.HO-Apoe.NC, Apoe.NC-Apoe.HO,...
                 Apoe.HE-Apoe.HO, Apoe.HO-Apoe.HE,...
                 (0.5*Apoe.HE+0.5*Apoe.NC) - Apoe.HO,...
                 Apoe.HO - (0.5*Apoe.HE+0.5*Apoe.NC),...
                 (0.5*Apoe.HE+0.5*Apoe.HO) - Apoe.NC,...
                 Apoe.NC - (0.5*Apoe.HE+0.5*Apoe.HO),...
                 DX.AD-DX.CN, DX.CN-DX.AD,...
                 DX.LMCI-DX.CN, DX.CN-DX.LMCI,...
                 DX.AD-DX.LMCI, DX.LMCI-DX.AD,...
                 DX.AD - (0.5*DX.LMCI + 0.5*DX.CN),...
                 DX.CN - (0.5*DX.AD + 0.5*DX.LMCI)
                 };


% List of files to save
save_file = {'-age', '+age', 'site', '-site'...
             '+male_-female', '+female_-male', '-volume', '+volume',...
             '+HE_-NC', '+NC_-HE', '+HO_-NC', '+NC_-HO', '+HO_-HE',...  
             '+HE_-HO', 'HE+NC-HO', '-HE-NC+HO', 'HE+HO-NC', '-HE-HO+NC',...
             '+AD-CN','-AD+CN','+MCI-CN','+CN-LMCI','+AD-MCI','-AD+MCI',...
             'CN+MCI-AD','CN-MCI+AD'};

% ADDITIVE MODEL
contrast_list_add = {apoe_int, -apoe_int, dx_int, -dx_int};
save_file_add = {'NC_HE_HO', '-NC_HE_HO', 'AD_MCI_CN', '-AD_MCI_CN'};
   
%% Left
         
% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo/');
mkdir(exp_dir);


X = []; % Part of the model we do not want to correct for
save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, apoe, corr)

% ADDITIVE MODEL
X = []; % Part of the model we do not want to correct for
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, age, apoe, corr_int)

%% Right

% List of contrasts

contrast_list = {-age, age, site, -site...
                 Gender.M-Gender.F, Gender.F-Gender.M,...
                 -volume_r, volume_r,...
                 Apoe.HE-Apoe.NC, Apoe.NC-Apoe.HE,...
                 Apoe.HO-Apoe.NC, Apoe.NC-Apoe.HO,...
                 Apoe.HE-Apoe.HO, Apoe.HO-Apoe.HE,...
                 (0.5*Apoe.HE+0.5*Apoe.NC) - Apoe.HO,...
                 Apoe.HO - (0.5*Apoe.HE+0.5*Apoe.NC),...
                 (0.5*Apoe.HE+0.5*Apoe.HO) - Apoe.NC,...
                 Apoe.NC - (0.5*Apoe.HE+0.5*Apoe.HO),...
                 DX.AD-DX.CN, DX.CN-DX.AD,...
                 DX.LMCI-DX.CN, DX.CN-DX.LMCI,...
                 DX.AD-DX.LMCI, DX.LMCI-DX.AD,...
                 DX.AD + (0.5*DX.LMCI - 0.5*DX.CN),...
                 DX.CN - (0.5*DX.AD + 0.5*DX.LMCI)
                 };


% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, apoe, corr)

% ADDITIVE MODEL
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, apoe, corr_int)
