%% Load OR ratings
% We assume that data is already loaded

apoe_allele = strcat('A',cellstr(string(data.APOE_allele)));
APOE_allele = term(apoe_allele);

apoe_OR = data.APOE_OR;
APOE_OR = term(apoe_OR);

apoe_OR_sq = data.APOE_OR .* data.APOE_OR;%  - mean(data.age .* data.age);
APOE_OR_sq = term(apoe_OR_sq);
