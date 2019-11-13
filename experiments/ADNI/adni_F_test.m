% Baseline experiments on rigid alineation
% Test left and right hippocampus separately
%% Do multivariate statistical surface analysis with freestat, akin to the paper
global MAXMEM; MAXMEM=4096;
clear all
close all
addpath('surfstat/')
addpath('utils/')

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
Data_loading_ADNI;

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

model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe + DX;
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe + DX;

slm_l = SurfStatLinMod(t_l, model_l, avg_l);
slm_r = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int + DX_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int + DX_int;

slm_l_add = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add = SurfStatLinMod(t_r, model_r_add, avg_r);

% Interaction model (Age x Apoe)
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe + DX + Age*Apoe + Agesq*Apoe;
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe + DX + Age*Apoe + Agesq*Apoe;

slm_l_ageapoe = SurfStatLinMod(t_l, model_l, avg_l);
slm_r_ageapoe = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int + DX_int + Age*Apoe_int + Agesq*Apoe_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int + DX_int + Age*Apoe_int + Agesq*Apoe_int;

slm_l_add_ageapoe = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add_ageapoe = SurfStatLinMod(t_r, model_r_add, avg_r);

% Interaction model (Age x DX)
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe + DX + Age*DX + Agesq*DX;
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe + DX + Age*DX + Agesq*DX;

slm_l_agedx = SurfStatLinMod(t_l, model_l, avg_l);
slm_r_agedx = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int + DX_int + Age*DX_int + Agesq*DX_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int + DX_int + Age*DX_int + Agesq*DX_int;

slm_l_add_agedx = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add_agedx = SurfStatLinMod(t_r, model_r_add, avg_r);

% Interaction model (Apoe x DX)
model_l = 1 + Age + Agesq + Gender + Volume_l + Apoe + DX + Apoe*DX;
model_r = 1 + Age + Agesq + Gender + Volume_r + Apoe + DX + Apoe*DX;

slm_l_apoedx = SurfStatLinMod(t_l, model_l, avg_l);
slm_r_apoedx = SurfStatLinMod(t_r, model_r, avg_r);

% ADDITIVE
model_r_add = 1 + Age + Agesq + Gender + Volume_r + Apoe_int + DX_int + Apoe_int*DX_int;
model_l_add = 1 + Age + Agesq + Gender + Volume_l + Apoe_int + DX_int + Apoe_int*DX_int;

slm_l_add_apoedx = SurfStatLinMod(t_l, model_l_add, avg_l);
slm_r_add_apoedx = SurfStatLinMod(t_r, model_r_add, avg_r);

%% F TESTING
% For each test, save to disk the original map and the FDR corrected map.
% I will end up writing a function for this... 


%% ageapoe
slm_ageapoe_r = SurfStatF( slm_r_ageapoe, slm_r );
SurfStatView( slm_ageapoe_r.t(1,:), avg_r, 'F statistic');
saveas(gcf,strcat(exp_dir_base, 'fstat_r_ageapoe.png'))
close

SurfStatView( SurfStatQ( slm_ageapoe_r, [] ), avg_r, 'F statistic, right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_r_ageapoe_fdr.png'))
close

slm_add_r_ageapoe = SurfStatF( slm_r_add_ageapoe, slm_r_add );
SurfStatView( slm_add_r_ageapoe.t(1,:), avg_r, 'F statistic, additive right hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_ageapoe.png'))
close

SurfStatView( SurfStatQ( slm_add_r_ageapoe, [] ), avg_r, 'F statistic, additive right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_ageapoe_fdr.png'))
close

slm_ageapoe_l = SurfStatF( slm_l_ageapoe, slm_l );
SurfStatView( slm_ageapoe_l.t(1,:), avg_l, 'F statistic, left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_ageapoe.png'))
close

SurfStatView( SurfStatQ( slm_ageapoe_l, [] ), avg_l, 'F statistic, left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_ageapoe_fdr.png'))
close

slm_add_l_ageapoe = SurfStatF( slm_l_add_ageapoe, slm_l_add );
SurfStatView( slm_add_l_ageapoe.t(1,:), avg_l, 'F statistic, additive left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_ageapoe.png'))
close

SurfStatView( SurfStatQ( slm_add_l_ageapoe, [] ), avg_l, 'F statistic, additive left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_ageapoe_fdr.png'))
close

%% apoedx
slm_apoedx_r = SurfStatF( slm_r_apoedx, slm_r );
SurfStatView( slm_apoedx_r.t(1,:), avg_r, 'F statistic');
saveas(gcf,strcat(exp_dir_base, 'fstat_r_apoedx.png'))
close

SurfStatView( SurfStatQ( slm_apoedx_r, [] ), avg_r, 'F statistic, right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_r_apoedx_fdr.png'))
close

slm_add_r_apoedx = SurfStatF( slm_r_add_apoedx, slm_r_add );
SurfStatView( slm_add_r_apoedx.t(1,:), avg_r, 'F statistic, additive right hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_apoedx.png'))
close

SurfStatView( SurfStatQ( slm_add_r_apoedx, [] ), avg_r, 'F statistic, additive right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_apoedx_fdr.png'))
close

slm_apoedx_l = SurfStatF( slm_l_apoedx, slm_l );
SurfStatView( slm_apoedx_l.t(1,:), avg_l, 'F statistic, left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_apoedx.png'))
close

SurfStatView( SurfStatQ( slm_apoedx_l, [] ), avg_l, 'F statistic, left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_apoedx_fdr.png'))
close

slm_add_l_apoedx = SurfStatF( slm_l_add_apoedx, slm_l_add );
SurfStatView( slm_add_l_apoedx.t(1,:), avg_l, 'F statistic, additive left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_apoedx.png'))
close

SurfStatView( SurfStatQ( slm_add_l_apoedx, [] ), avg_l, 'F statistic, additive left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_apoedx_fdr.png'))
close

%% agedx
slm_agedx_r = SurfStatF( slm_r_agedx, slm_r );
SurfStatView( slm_agedx_r.t(1,:), avg_r, 'F statistic');
saveas(gcf,strcat(exp_dir_base, 'fstat_r_agedx.png'))
close

SurfStatView( SurfStatQ( slm_agedx_r, [] ), avg_r, 'F statistic, right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_r_agedx_fdr.png'))
close

slm_add_r_agedx = SurfStatF( slm_r_add_agedx, slm_r_add );
SurfStatView( slm_add_r_agedx.t(1,:), avg_r, 'F statistic, additive right hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_agedx.png'))
close

SurfStatView( SurfStatQ( slm_add_r_agedx, [] ), avg_r, 'F statistic, additive right hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_r_agedx_fdr.png'))
close

slm_agedx_l = SurfStatF( slm_l_agedx, slm_l );
SurfStatView( slm_agedx_l.t(1,:), avg_l, 'F statistic, left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_agedx.png'))
close

SurfStatView( SurfStatQ( slm_agedx_l, [] ), avg_l, 'F statistic, left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_l_agedx_fdr.png'))
close

slm_add_l_agedx = SurfStatF( slm_l_add_agedx, slm_l_add );
SurfStatView( slm_add_l_agedx.t(1,:), avg_l, 'F statistic, additive left hipp' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_agedx.png'))
close

SurfStatView( SurfStatQ( slm_add_l_agedx, [] ), avg_l, 'F statistic, additive left hipp, FDR' );
saveas(gcf,strcat(exp_dir_base, 'fstat_add_l_agedx_fdr.png'))
close

