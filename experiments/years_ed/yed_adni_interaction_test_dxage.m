% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;
clear all
close all
addpath('../surfstat/')
addpath('../utils/')
set(groot,'defaultFigureVisible','off')

%% 1. Load the data paths and covariates
meshes_l = "";
meshes_r = "";
csv = "";

avg_l_path = "";
avg_r_path = "";
avg_l_vtk = "";
avg_r_vtk = "";

exp_dir_base = "";
mkdir(exp_dir_base);

data = readtable(csv);
% Number of subjects
N = size(data,1);

%% Load covariates
Data_loading_ADNI;

%% Models
% Left hippocampus

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

model_l = 1 + Age + Agesq + Gender + Site + DX + Apoe + Yed + Age*DX + Agesq*DX;
model_r = 1 + Age + Agesq + Gender + Site + DX + Apoe + Yed + Age*DX + Agesq*DX;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Site + DX_int + Yed + Apoe_int + Age*DX_int + Agesq*DX_int;
model_l_add = 1 + Age + Agesq + Gender + Site + DX_int + Yed + Apoe_int + Age*DX_int + Agesq*DX_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);

% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and 
% the three diagnosis.

% for the additive part, is different, as we only have a covariate.

% CONFIRMAR
corr = [1 4 5 6 10 11 12 13];
corr_int = [1 4 5 6 8 9];


%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)

%% Contrasts

% NO INTERACTIONS
% List of contrasts
contrast_list = {-age, age, -agesq, agesq, site, -site, yed, -yed,...
                 Gender.Male-Gender.Female, Gender.Female-Gender.Male,...%-volume_l, volume_l,...
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
save_file = {'-age', '+age', '-agesq', 'agesq', 'site', '-site', 'yed', '-yed',...
             '+male_-female', '+female_-male',...% '-volume', '+volume',...
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

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, dx, corr)

% ADDITIVE MODEL
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, age, dx, corr_int)


%% INTERACTIONS TIME

contrast_list_int = {age.*DX.LMCI-age.*DX.CN, age.*DX.CN-age.*DX.LMCI,...
                 age.*DX.AD-age.*DX.CN, age.*DX.CN-age.*DX.AD,...
                 age.*DX.LMCI-age.*DX.AD, age.*DX.AD-age.*DX.LMCI,...
                 age.*(0.5*DX.LMCI+0.5*DX.CN) - age.*DX.AD,...
                 age.*DX.AD - age.*(0.5*DX.LMCI+0.5*DX.CN),...
                 age.*(0.5*DX.LMCI+0.5*DX.AD) - age.*DX.CN,...
                 age.*DX.CN - age.*(0.5*DX.LMCI+0.5*DX.AD),...
                 agesq.*DX.LMCI-agesq.*DX.CN, agesq.*DX.CN-agesq.*DX.LMCI,...
                 agesq.*DX.AD-agesq.*DX.CN, agesq.*DX.CN-agesq.*DX.AD,...
                 agesq.*DX.LMCI-agesq.*DX.AD, agesq.*DX.AD-agesq.*DX.LMCI,...
                 agesq.*(0.5*DX.LMCI+0.5*DX.CN) - agesq.*DX.AD,...
                 agesq.*DX.AD - agesq.*(0.5*DX.LMCI+0.5*DX.CN),...
                 agesq.*(0.5*DX.LMCI+0.5*DX.AD) - agesq.*DX.CN,...
                 agesq.*DX.CN - agesq.*(0.5*DX.LMCI+0.5*DX.AD)};
    
% List of files to save
% 
save_file_int = {'age+MCI_-CN', 'age+CN_-MCI', 'age+AD_-CN', 'age+CN_-AD', 'age+AD_-MCI',...  
             'age+MCI_-AD', 'ageMCI+CN-AD', 'age-MCI-CN+AD', 'ageMCI+AD-CN', 'age-MCI-AD+CN',...
             'agesq+MCI_-CN', 'agesq+CN_-MCI', 'agesq+AD_-CN', 'agesq+CN_-AD', 'agesq+AD_-MCI',...  
             'agesq+MCI_-AD', 'agesqMCI+CN-AD', 'agesq-MCI-CN+AD', 'agesqMCI+AD-CN', 'agesq-MCI-AD+CN'};

% ADDITIVE
contrast_list_int_add = {age.*dx_int, -age.*dx_int, agesq.*dx_int, -agesq.*dx_int};
save_file_int_add = {'ageCN_MCI_AD', 'age-CN_MCI_AD', 'agesqCN_MCI_AD', 'agesq-CN_MCI_AD'};


% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo_ageapoe/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list_int,save_file_int,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, dx, corr)

% ADDITIVE
save_exp_to_disk(contrast_list_int_add,save_file_int_add,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, age, dx, corr_int)

%% Right

contrast_list = {-age, age, -agesq, agesq, site, -site, yed, -yed,...
                 Gender.Male-Gender.Female, Gender.Female-Gender.Male,...% -volume_r, volume_r,...
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

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, dx, corr)

% ADDITIVE MODEL
save_exp_to_disk(contrast_list_add,save_file_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, dx, corr_int)


% Interactions time

exp_dir = strcat(exp_dir_base, 'righthippo_ageapoe/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list_int,save_file_int,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, dx, corr)

% ADDITIVE
save_exp_to_disk(contrast_list_int_add,save_file_int_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, dx, corr_int)

