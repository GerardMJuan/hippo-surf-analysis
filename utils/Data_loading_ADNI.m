
% Age and age squared, not centered at mean
age = data.AGE - mean(data.AGE);
Age = term(age);

agesq = data.AGE .* data.AGE - mean(data.AGE .* data.AGE);
Agesq = term(agesq);

% Gender, no further things
gender = data.PTGENDER;
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

apoe_intsq = data.apoe_int .* data.apoe_int;
Apoe_intsq = term(apoe_intsq);

% Volume
volume_l = data.vol_l;
Volume_l = term(volume_l);

volume_r = data.vol_r;
Volume_r = term(volume_r);

% DX
dx = data.DX_bl;
DX = term(dx);

dx_int = data.dx_int;
DX_int = term(dx_int);

% Years of education
yed = data.PTEDUCAT;
Yed = term(yed);

yedsq = data.PTEDUCAT .* data.PTEDUCAT;
Yedsq = term(yedsq);

% Load the images
id = data.PTID;
subdir_l = strcat(meshes_l, string(id),'_l.obj');
subdir_r = strcat(meshes_r, string(id),'_r.obj');

%% Put the meshes in a format that is legible by surfstat

mesh_l = SurfStatReadSurf(subdir_l);
mesh_r = SurfStatReadSurf(subdir_r);

avg_l = SurfStatReadSurf([avg_l_path]);
avg_r = SurfStatReadSurf([avg_r_path]);

mask_l = true(1, size(avg_l.coord, 2));
mask_r = true(1, size(avg_r.coord, 2));
