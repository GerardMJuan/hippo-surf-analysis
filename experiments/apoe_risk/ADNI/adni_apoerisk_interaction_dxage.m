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

avg_l_path = '';
avg_r_path = '';
avg_l_vtk = '';
avg_r_vtk = '';

exp_dir_base = '';
mkdir(exp_dir_base);

data = readtable(csv);
% Number of subjects
N = size(data,1);

%% Load covariates
Data_loading_ADNI;
% Load the apoe OR
or_loading;

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

model_l = 1 + Age + Agesq + Gender + Site + DX + APOE_OR + Age*DX + Agesq*DX;
model_r = 1 + Age + Agesq + Gender + Site + DX + APOE_OR + Age*DX + Agesq*DX;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);


% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Site + DX_int + APOE_OR + Age*DX_int + Agesq*DX_int;
model_l_add = 1 + Age + Agesq + Gender + Site + DX_int + APOE_OR + Age*DX_int + Agesq*DX_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);


% for the additive part, is different, as we only have a covariate.

% CONFIRMAR
corr = [1 4 5 9 10 11];
corr_int = [1 4 5 7];
corr = [];
corr_int = [];

%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)

%% Contrasts

% NO INTERACTIONS
% List of contrasts
contrast_list = {-age, age, -agesq, agesq, site, -site,...
                 Gender.Male-Gender.Female, Gender.Female-Gender.Male,...%-volume_l, volume_l,...
                 apoe_OR, -apoe_OR,...
                 DX.AD-DX.CN, DX.CN-DX.AD,...
                 DX.LMCI-DX.CN, DX.CN-DX.LMCI,...
                 DX.AD-DX.LMCI, DX.LMCI-DX.AD,...
                 DX.AD - (0.5*DX.LMCI + 0.5*DX.CN),...
                 DX.CN - (0.5*DX.AD + 0.5*DX.LMCI)
                 };

% List of files to save
save_file = {'-age', '+age', '-agesq', 'agesq', 'site', '-site'...
             '+male_-female', '+female_-male',...% '-volume', '+volume',...
             'apoe_or', '-apoe_or',...
             '+AD-CN','-AD+CN','+MCI-CN','+CN-LMCI','+AD-MCI','-AD+MCI',...
             'CN+MCI-AD','CN-MCI+AD'};


%% Left
% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo/');
mkdir(exp_dir);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, dx, corr)

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

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, dx, corr)


% Interactions time

exp_dir = strcat(exp_dir_base, 'righthippo_ageapoe/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list_int,save_file_int,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, dx, corr)

% ADDITIVE
save_exp_to_disk(contrast_list_int_add,save_file_int_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, dx, corr_int)

