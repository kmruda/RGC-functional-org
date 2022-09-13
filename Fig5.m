
%% Get temporal RFs in frequency space and extract peak frequency

% all_raw_tcs: structure of temporal RFs for each cell type/light level

types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
ndf_names = {'NDF5';'NDF4';'NDF3';'NDF2';'NDF1';'NDF0'};

% normalize the raw time courses
for type = 1:length(types_interest)
    for j=1:length(all_ndf_names)
        these_tcs = all_raw_tcs.(types_interest{type}).(all_ndf_names{j});
        for rgc = 1:size(these_tcs,2)
            this_norm_tc = these_tcs(:,rgc)/norm(these_tcs(:,rgc));
            all_normd_tcs.(types_interest{type}).(all_ndf_names{j})(:,rgc) = this_norm_tc;
        end
    end
end

% make time vectors 
intervals = [4 2 1 1 1 1];
nkts = [30 30 30 30 30 30];
time_intervals=intervals/60.35;

% get power spectra
fs = 60./intervals;
for type=1:length(types_interest)
    for j=1:length(all_ndf_names)
        these_tcs = all_raw_tcs.(types_interest{type}).(all_ndf_names{j});
        samplerate = fs(j);
        for rgc=1:size(these_tcs,2)
            temp_tc = all_normd_tcs.(types_interest{type}).(all_ndf_names{j})(:,rgc);
            [powerspec_xvalues, mean_powerspec, fft_signal] = PowerSpectrumFinder(temp_tc',samplerate) ;
            freq_ranges.(all_ndf_names{j}) = powerspec_xvalues;
            norm_ps = mean_powerspec/max(mean_powerspec);
            PS_norm.(types_interest{type}).(all_ndf_names{j})(:,rgc) = norm_ps;
        end
    end
end

% get peak frequency
for type=1:length(types_interest)
    for j=1:length(all_ndf_names)
        for rgc=1:size(PS_norm.(types_interest{type}).(all_ndf_names{j}),2)
            if ~isnan(PS_norm.(types_interest{type}).(all_ndf_names{j})(1,rgc))
                [maxs,inds] = max(PS_norm.(types_interest{type}).(all_ndf_names{j})(:,rgc));
                pk_freqs.(types_interest{type}).(all_ndf_names{j})(rgc) = freq_ranges.(all_ndf_names{j})(inds);
            else
                pk_freqs.(types_interest{type}).(all_ndf_names{j})(rgc) = nan;
            end
        end
        avg_pk_freqs(type,j) = nanmean(pk_freqs.(types_interest{type}).(all_ndf_names{j}));
        std_pk_freqs(type,j) = nanstd(pk_freqs.(types_interest{type}).(all_ndf_names{j}));
        num_pk_freqs(type,j) = sum(~isnan(pk_freqs.(types_interest{type}).(all_ndf_names{j})));
    end
end