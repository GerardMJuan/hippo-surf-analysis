% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
% This is new!
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

%%% 

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

% ADDITIVE
% BUT USING OR
model_r_add = 1 + Age + Gender + APOE_OR + Yed;
model_l_add = 1 + Age + Gender + APOE_OR + Yed;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);


% Corr is a matrix with the indexes of the design matrix that need to be 
% corrected for the age/apoe interaction figure
% In this case, it represents the intercept, gender (M and F) and volume.
% if it contains nothing, then it interprets that there is no interaction
% and do not compute the figure (see utils/save_peaks_figure.m)
corr = [];


%% Save uncorrected tests to disk
% Save, for each contrast:
% The uncorrected image
% The corrected image
% the .vtk with the uncorrected maps on it, somehow (call to python?)
% List of contrasts

contrast_list = {-age, age,...% volume_l, -volume_l,...
                 Gender.M-Gender.F, Gender.F-Gender.M, yed, -yed,...
                 apoe_OR, -apoe_OR
                 };

% List of files to save
%
save_file = {'-age', '+age',...% 'volume', '-volume'...
             '+male_-female', '+female_-male','yed', '-yed',...
             'apoe_or', '-apoe_or'};

% output_dir
exp_dir = strcat(exp_dir_base, 'lefthippo/');
mkdir(exp_dir);

% ADDITIVE MODEL
save_exp_to_disk(contrast_list,save_file,exp_dir,avg_l,slm_l_add,avg_l_vtk, t_l, age, apoe, corr)

%% Right

% output_dir
exp_dir = strcat(exp_dir_base, 'righthippo/');
mkdir(exp_dir);

% ADDITIVE MODEL
save_exp_to_disk(contrast_list,save_file,exp_dir,avg_r,slm_r_add,avg_r_vtk, t_r, age, apoe, corr)
