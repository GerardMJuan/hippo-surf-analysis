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


% ADDITIVE
model_r_add = 1 + Age + Gender + Site + APOE_OR + DX_int +APOE_OR*DX_int + APOE_OR_sq*DX_int;
model_l_add = 1 + Age + Gender + Site + APOE_OR + DX_int +APOE_OR*DX_int + APOE_OR_sq*DX_int;
slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);

% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and 
% the three diagnosis.

% for the additive part, is different, as we only have a covariate.

% CONFIRMAR
corr = [];
corr_int = [1 2 3 4];


%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)

%% Left


% INTERACTIONS TIME
% ONLY ADDITIVE
% ADDITIVE
contrast_list_int_add = {apoe_OR.*dx_int, -apoe_OR.*dx_int, apoe_OR_sq.*dx_int, -apoe_OR_sq.*dx_int};
save_file_int_add = {'apoeDX', '-apoeDX', 'apoesqDX', '-apoesqDX'};

% HACK FINS QUE HAGI CORREGIT EL GRAFIC PER AQUEST CAS
corr_int = [];
% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo_ageapoe/');
mkdir(exp_dir);

% % Call save function
save_exp_to_disk(contrast_list_int_add,save_file_int_add,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, apoe, dx, corr_int)

%% Right

% INTERACTIONS TIME

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo_ageapoe/');
mkdir(exp_dir);

% % Call save function
save_exp_to_disk(contrast_list_int_add,save_file_int_add,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, apoe, dx, corr_int)
