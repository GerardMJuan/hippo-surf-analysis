%% Load covariates
% No hi ha raó per la qual s'ha de dividir en 2 cada cosa, tenint en compte
% que tenen la mateixa llargada i a més contenen la mateixa informació
% Age
% Tenia una raó però no la recordo. Era bona? No sé...
age = data.age - mean(data.age);
Age = term(age);

agesq = data.age .* data.age - mean(data.age .* data.age);
Agesq = term(agesq);

% Gender, no further
gender = data.gender;
Gender = term(gender);

% site
site = data.SITE;
Site = term(site);

% Apoe, genotype model
% Divide into two dummy regressors
apoe = data.apoe_cat;
Apoe = term(apoe);

apoe_int = data.apoe_int;
Apoe_int = term(apoe_int);

% DX
dx = data.DX_bl;
DX = term(dx);

dx_int = data.dx_int;
DX_int = term(dx_int);

% Volume
% tot i que no el farem servir
volume_l = data.vol_l;
Volume_l = term(volume_l);

volume_r = data.vol_r;
Volume_r = term(volume_r);

% Years of education
yed = data.education_years;
Yed = term(yed);

yedsq = data.education_years .* data.education_years;
Yedsq = term(yedsq);

% Load the images
id = data.PID;
subdir_l = strcat(meshes_l,'out_mesh_', string(id),'_l.obj');
subdir_r = strcat(meshes_r,'out_mesh_', string(id),'_r.obj');

%% Put the meshes in a format that is legible by surfstat

mesh_l = SurfStatReadSurf(subdir_l);
mesh_r = SurfStatReadSurf(subdir_r);

avg_l = SurfStatReadSurf([avg_l_path]);
avg_r = SurfStatReadSurf([avg_r_path]);

mask_l = true(1, size(avg_l.coord, 2));
mask_r = true(1, size(avg_r.coord, 2));
