
%% Get nonlinear index from fit SNLs

% fit_params: parameter values for fitted SNLs from datarun structure

types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
ndf_names = {'NDF5';'NDF4';'NDF3';'NDF2';'NDF1';'NDF0'};

snl_range = 1.5;
num_points = 100;

for type = 1:length(types_interest)
    for j = 1:length(ndf_names)
        temp_indices = get_cell_indices(dataruns.(ndf_names{j}), IDZ_run.(types_interest{type}){j});
        for rgc = 1:length(temp_indices)
            fit_params = dataruns.(ndf_names{j}).stas.snls{temp_indices(rgc)}.fit_params;
            fitx = linspace(-1*snl_range, snl_range, num_points);
            fitfcn = fit_params.a * normcdf(fit_params.b * fitx - fit_params.c, 0, 1);        
            slope_at_zero = abs(diff(fitfcn(50:51))) / ((snl_range*2)./num_points);
            slope_at_max = max(diff(fitfcn));
            
            temp_nli = log(slope_at_max / slope_at_zero);

            dataruns.(ndf_names{j}).stas.snls{temp_indices(rgc)}.nli = temp_nli;
            if isempty(temp_nli)
                NL_index.(types_interest{type}){j}(rgc) = nan;
            else
                NL_index.(types_interest{type}){j}(rgc) = temp_nli;
            end
        end
    end
end


