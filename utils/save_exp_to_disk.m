function save_exp_to_disk(contrast_list, save_file, exp_dir, avg, slm_b, avg_vtk, Y, age, apoe, corr)
% Auxiliar function to not repeat code. Saves all the experiment contained
% in contrast_list to exp_dir.
% It also calls on save_peaks_figure, whenever we are testing for
% interactions and survive the correction
% y, age and apoe are used to do the figures
% corr is for correcting the input for the actual graphic
% In the end, save a csv with the changes 
    for i = 1:length(contrast_list)
        disp(save_file{i});
        % Uncorrected contrast
        slm = SurfStatT( slm_b, contrast_list{i});
        % Corrected contrast
        [ pval, peak, clus, clusid ] = SurfStatP( slm, [], 0.001 );

        % Qvalues
        qval = SurfStatQ( slm ); 
        SurfStatView( qval, avg, save_file{i} );
        saveas(gcf,strcat(exp_dir, save_file{i},'_fdr_corr.png'))
        close

        % normal vector for each point
        N = patchnormals(avg);
        
        % TODO: save also clus and peak to a csv o similars
        if ~isempty(clus)
            clus_out = strcat(exp_dir, save_file{i},'_clusters.txt');
            writetable(struct2table(clus), clus_out)
        end

        if ~isempty(peak)
            peak_out = strcat(exp_dir, save_file{i},'_peak.txt');
            writetable(struct2table(peak), peak_out)
        end
        
        %if we have an X to correct, means that there is a 
        if ~isempty(peak) && ~isempty(clus) && ~isempty(corr)
            % Call the function save_peaks_figure
            mkdir(strcat(exp_dir, save_file{i},'/'))
            %save_peaks_figure(Y, N, slm, pval, clus, peak,...
            %    clusid, age, apoe, strcat(exp_dir, save_file{i},'/'), corr)
        end
        
        % Visualize uncorrected and save to disk
        SurfStatView( slm.t, avg, save_file{i} );
        saveas(gcf,strcat(exp_dir, save_file{i},'_uncorr.png'))
        close

        % Save values to disk for them to be used in vtk to create a vtk later
        % We need to save:
        % 1: T-statistic (Hotelling's T)
        % The t statistic will need to have a different sign depending on
        % the orientation of the vector w.r.t. the normal of the point
        % 2: T values for x coordinate
        % 3: T values for y coordinate
        % 4: T values for z coordinate
                
        % Abans h otenia multiplicat per minus, why
        tstat = slm.t .* sign(dot(N',slm.ef./slm.sd));
        angle = atan2(norm(cross(N',slm.ef./slm.sd)),dot(N',slm.ef./slm.sd));
        
        csv_out = strcat(exp_dir, save_file{i},'_values.csv');
        to_csv_mat = [tstat', angle', slm.ef'];
        csvwrite(csv_out, to_csv_mat);

        % Call the function from python to create the volume
        vtk_out = strcat(exp_dir, save_file{i},'_Tvalues.vtk');
        system("python utils/create_vtk_from_matlab.py " + avg_vtk + " " +...
               csv_out + " " + vtk_out);

        % TEMPORAL
        csv_out = strcat(exp_dir, 'normals.csv');
        to_csv_mat = [tstat', angle', N];
        csvwrite(csv_out, to_csv_mat);
        vtk_out = strcat(exp_dir, 'normals.vtk');
        system("python utils/create_vtk_from_matlab.py " + avg_vtk + " " +...
               csv_out + " " + vtk_out);

           
        % Visualize corrected and save to disk
        SurfStatView( pval, avg, save_file{i} );
        saveas(gcf,strcat(exp_dir, save_file{i},'_corr.png'))
        close

        % Visualize and save cluster position
        SurfStatView( clusid, avg, save_file{i} );
        saveas(gcf,strcat(exp_dir, save_file{i},'_clusid.png'))
        close
    end
    
end