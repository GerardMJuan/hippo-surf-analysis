function save_peaks_figure(Y, N, slm, pval, clus, peak, clusid, age, apoe, out_dir, corr)
%SAVE_PEAKS_FIGURE Auxiliar function to create figures for significant
%peaks and clusters of a test. 
%   For a given interaction with APOE, we correct the results on 
% the surface by removing the effect of the covariates, create a csv with
% the relevant information, and run a python script that creates the
% corresponding figures
% Y is the response variable (distance to mean shape, normally)
% N are the normals for each vertex in order to normalize the direction
% pval contains the pvaules of the clusters and regions
% slm is the fitted model with all the covariates
% clus is the list of clusters found, obtained in SurfStatP
% Same for peaks and clusid
% age and apoe, information about those to make the graphics
% out_dir is the output directory 
% corr is list of inputs that need to be corrected
%
% Yadj = zeros(size(Y)); % Where the adjusted values will be saved
% % Correct response variable
% for i= 1:size(Y,3)
%     for j = 1:size(Y,2)
%         % Do the adjustment similar to the one did in SurfstatPlot
%         y = double(Y(:,j,i));
%         Yhat=slm.X*squeeze(slm.coef(:,j,i));
%         XGmat=double(X);
%         Yhatadj=XGmat*pinv(XGmat)*Yhat;
%         c=mean(y)-mean(Yhatadj);
%         Yhatadj=Yhatadj+c;
%         Yadj(:,j,i) = Yhatadj+y-Yhat;
%     end
% end

% Save info to disk
% We should save:
% 1. residuals of all interest zones
% 2. terms of the interaction (linear and no-linear, to draw the curve)
% use clusid to select only the vertices of importance and save it all,
% together with the other informarion

% For the mean of all the vertexs of the indicated cluster, and for the 

%% Save info to disk
% We should save:
% 1. residuals of all interest zones
% 2. terms of the interaction (linear and no-linear, to draw the curve)
% use clusid to select only the vertices of importance and save it all,
% together with the other informarion
for i = 1:length(clus.clusid)
    % For each cluster, compute:
    
    % the mean effect over all the pixels of the
    % cluster
    save_and_call(Y, clus.P(i), clusid == i, N, slm, strcat('mean_eff_',...
        int2str(i)), age, apoe, out_dir, corr);
    
    % save the peak
    save_and_call(Y, mean(squeeze(pval.P(1, peak.vertid(i)))), peak.vertid(i), N, slm, strcat('peak_',...
        int2str(i)), age, apoe, out_dir, corr);

end
end

function save_and_call(Yadj, pval, index, N, slm, intfix, age, apoe, out_dir, corr)
%SAVE_AND_CALL Auxiliar function to not repeat code. 
%   Given a list of vertices and an output of a SurfStat correction,
%   outputs a csv with all the information and runs the python code
% Yadj: adjusted response variables
% index: indexs of the parts of the cluster that are relevant.
% N: Normals of the points on the hipppocampus surface to capture the area
% slm: SurfStat slm, output of a SurfstatP.
% age and apoe, information about those to make the graphics
% intfix: intfix name put to 
% corr is the list of inputs that need to be corrected
    % Get normals of the index
    N_ind = N(index,:); 

    Yhat = zeros(size(Yadj)); % Where the adjusted values will be saved
    for i= 1:size(Yhat,3)
        Yhat(:,index,i) = Yadj(:,index,i) - slm.X(:, corr)*squeeze(slm.coef(corr,index,i));
    end
    % Save X, Y and Z, in parallel
    % Find the mean and the norm

    % Agafo la mitjana sobre tota l'area de els indexs passats
    vert = Yhat(:,index, :);
    vert_mean = squeeze(mean(vert, 2));
    N_mean = squeeze(mean(N_ind, 1));

    vert_x = vert_mean(:,1);
    vert_y = vert_mean(:,2);
    vert_z = vert_mean(:,3);

    N_mean_repmat = repmat(N_mean, [length(vert_mean) 1]);

    % La norm de cada un, amb la mean 
    vert_norm = vecnorm(vert_mean, 2, 2) .* -sign(dot(vert_mean',N_mean_repmat'))';
    
    % The linear and quadratic effect of the age, save it in the csv too
    % and the pvalue of the cluster
    % THIS IS HARDCODED, SHOULD BE CHANGED
    intercept = vecnorm(squeeze(mean(slm.coef(1, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(1, index, :), 2)),N_mean));
    agelin = vecnorm(squeeze(mean(slm.coef(2, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(2, index, :), 2)),N_mean));
    agequad = vecnorm(squeeze(mean(slm.coef(3, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(3, index, :), 2)),N_mean));
    apoecof = vecnorm(squeeze(mean(slm.coef(7, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(7, index, :), 2)),N_mean));
    lin = vecnorm(squeeze(mean(slm.coef(8, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(8, index, :), 2)),N_mean));
    quad = vecnorm(squeeze(mean(slm.coef(9, index, :))), 2) * -sign(dot(squeeze(mean(slm.coef(9, index, :), 2)),N_mean));
    % pval = mean(squeeze(pval.P(1, index)));
    
    % values = [0,0,0,0,0,0,pval];
    
    values = [intercept, agelin, agequad, apoecof, lin, quad, pval];
    values = [values, zeros(1, length(age) - length(values))]';

    csv_out = strcat(out_dir,intfix,'.csv');
    writetable(cell2table([num2cell(vert_norm), num2cell(age), apoe, num2cell(values)]),csv_out,'writevariablenames',0);
    
    % Save also x, y, z
    %csv_out_x = strcat(out_dir,intfix,'_x.csv');
    %writetable(cell2table([num2cell(vert_x), num2cell(age), apoe, num2cell(values)]),csv_out_x,'writevariablenames',0);
    %csv_out_y = strcat(out_dir,intfix,'_y.csv');
    %writetable(cell2table([num2cell(vert_y), num2cell(age), apoe, num2cell(values)]),csv_out_y,'writevariablenames',0);
    %csv_out_z = strcat(out_dir,intfix,'_z.csv');
    %writetable(cell2table([num2cell(vert_z), num2cell(age), apoe, num2cell(values)]),csv_out_z,'writevariablenames',0);

    fig_out = strcat(out_dir,intfix,'.png');
    system("python utils/create_interaction_figure.py " + csv_out + " " +...
           fig_out);  
    % Create also figures for x y z
   
    %fig_out = strcat(out_dir,intfix,'_x.png');
    %system("python utils/create_interaction_figure.py " + csv_out_x + " " +...
    %       fig_out);  
    %fig_out = strcat(out_dir,intfix,'_y.png');
    %system("python utils/create_interaction_figure.py " + csv_out_y + " " +...
    %       fig_out);  
    %fig_out = strcat(out_dir,intfix,'_z.png');
    %system("python utils/create_interaction_figure.py " + csv_out_z + " " +...
    %       fig_out);  

end

