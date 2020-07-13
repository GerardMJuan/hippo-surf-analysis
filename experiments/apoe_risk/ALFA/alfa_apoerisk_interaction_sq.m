% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;
clear all
close all
addpath('surfstat/')
addpath('utils/')
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

%% Load the data
Data_loading;
% Load the apoe OR
or_loading;

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

model_r = 1 + Age + Agesq + Gender + Yed + APOE_OR + Age*APOE_OR + Agesq*APOE_OR;
model_l = 1 + Age + Agesq + Gender + Yed + APOE_OR + Age*APOE_OR + Agesq*APOE_OR;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);


% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and volume.
corr = [1 4 5 6];

%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)

%% Left

contrast_list = {-age, age, yed, -yed,...% volume_l, -volume_l,...
                 Gender.M-Gender.F, Gender.F-Gender.M,...
                 apoe_OR, -apoe_OR};
    
% List of files to save
save_file = {'-age', '+age', 'yed', '-yed',...% '+volume', '-volume',...
             '+male_-female', '+female_-male',... % '-yed', 'yed',...
             'APOE_OR', 'APOE_OR'};


% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo/');
mkdir(exp_dir);

save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, apoe_allele, corr)



% ADDITIVE
contrast_list_int = {age.*apoe_OR, -age.*apoe_OR, agesq.*apoe_OR, -agesq.*apoe_OR};
save_file_int = {'ageNC_HE_HO', 'age-NC_HE_HO', 'agesqNC_HE_HO', 'agesq-NC_HE_HO'};

 
% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo_ageapoe/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list_int,save_file_int,exp_dir,avg_l,slm_l,avg_l_vtk, t_l, age, apoe_allele, corr)

%% Right

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, apoe_allele, corr)


%% INTERACTIONS TIME

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo_ageapoe/');
mkdir(exp_dir);

% Call save function
save_exp_to_disk(contrast_list_int,save_file_int,exp_dir,avg_r,slm_r,avg_r_vtk, t_r, age, apoe_allele, corr)