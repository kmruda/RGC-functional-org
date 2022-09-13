
%% Extract parameters from temporal RFs

% all_raw_tcs: structure of temporal RFs for each cell type/light level

types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
all_ndf_names = {'NDF5';'NDF4';'NDF3';'NDF2';'NDF1';'NDF0'};

% normalize the raw time courses
for type = 1:length(types_interest)
    for j=1:length(all_ndf_names)
        these_tcs = all_raw_tcs.(types_interest{type}).(all_ndf_names{j});
        for rgc = 1:size(these_tcs,2)
            factr = sign(type-3.5); 
            this_tc = factr*these_tcs(:,rgc);
            this_norm_tc = factr*these_tcs(:,rgc)/max(this_tc);
            all_normd_tcs.(types_interest{type}).(all_ndf_names{j})(:,rgc) = this_norm_tc;
        end
    end
end

% make time vectors
intervals = [4 2 1 1 1 1];
nkts = [30 30 30 30 30 30];
time_intervals=intervals/60.35;
for j=1:length(all_ndf_names)
    time_vectors{j}=-(nkts(j)-1)*time_intervals(j):time_intervals(j):0;       
end

% make them all ONs
for type = 1:3
    for j=1:length(all_ndf_names)
        all_raw_tcs.(types_interest{type}).(all_ndf_names{j}) = -1*all_raw_tcs.(types_interest{type}).(all_ndf_names{j});
    end
end

% subsample and extract parameters
lnkts = [291 291 291 291 291 291];
for type = 1:length(types_interest)
    for j = 1:length(all_ndf_names)
        finer_time_vectors{j} = -(nkts(j)-1)*time_intervals(j):time_intervals(j)/10:0;        
        for rgc = 1:size(all_raw_tcs.(types_interest{type}).(all_ndf_names{j}),2) 
            if ~isnan(all_raw_tcs.(types_interest{type}).(all_ndf_names{j})(1,rgc))
                fine_tc = interp1(time_vectors{j},all_raw_tcs.(types_interest{type}).(all_ndf_names{j})(:,rgc),finer_time_vectors{j});
                norm_fine_tc = factr*fine_tc/max(fine_tc); 
                [TCParams] = time_course_parameters(norm_fine_tc', 0, finer_time_vectors{j});
                tc_params.(types_interest{type}).(all_ndf_names{j}){rgc} = TCParams;
            end
        end
    end
end

% reorganize
metric_names={'stddev';'var';'ave';'minval';'mintim';'maxval';'maxtim';'norm';'areaabs';'firstpeaktime';'secondpeaktime';...
    'zerocrossing';'pk_dot'};
for type=1:length(types_interest)
    for j=1:length(all_ndf_names)
        for tp=1:length(metric_names)
            for rgc = 1:size(all_raw_tcs.(types_interest{type}).(all_ndf_names{j}),2) 
                if ~isnan(all_raw_tcs.(types_interest{type}).(all_ndf_names{j})(1,rgc))
                    this_param.(types_interest{type}){tp}(j,rgc) = tc_params.(types_interest{type}).(all_ndf_names{j}){rgc}.(metric_names{tp});
                else
                    this_param.(types_interest{type}){tp}(j,rgc) = nan;
                end
                
            end
        end
    end
end
for type=1:length(types_interest)
    for tp=1:length(metric_names)
        this_param.(types_interest{type}){tp}(this_param.(types_interest{type}){tp}==0)=nan;
    end
end

%%
%Calculation of Parameters from Temporal Receptive Field
function[TCParams] = time_course_parameters(timecourse, timevector)
% from Sneha Ravi
% Input: Time course matrix (rows - Tc values, columns - each cell) 
%        If correction is needed because of noisy time courses (1 or 0)
% Returns structure with parameters calculated from each cell's time course in this order:
% 1. Standard Deviation
% 2. Variance
% 3. Mean
% 4. Minimum Value
% 5. Minimum Value Time
% 6. Peak Value
% 7. Peak Value Time
% 8. Norm
% 9. Absolute Area Under Curve
% 10. Time to first peak
% 11. Time to second peak
% 12.Time of zero crossing
% 13. Total Area Under Curve / Absolute Area Under Curve (Degree of Transience)


TCParams = struct;

TCParams.stddev = std(timecourse);
TCParams.var = var(timecourse);
TCParams.ave = mean(timecourse);
[TCParams.minval TCParams.mintim] =  min(timecourse);
[TCParams.maxval TCParams.maxtim] =  max(timecourse);

time_interval = mean(diff(timevector));
for i = 1:size(timecourse, 2) 
    TCParams.norm(1, i) = norm( timecourse(:,i)); 
    TCParams.areaabs(1,i) = trapz(abs(timecourse(:,i)))*time_interval; 
    TCParams.area(1,i) = trapz(timecourse(:,i))*time_interval;
    if(TCParams.mintim(1,i) < TCParams.maxtim(1,i)) %if on 
        TCParams.firstpeaktime(1,i) = TCParams.maxtim(1,i);
        TCParams.secondpeaktime(1,i)= TCParams.mintim(1,i); 
    else %if off
        TCParams.firstpeaktime(1,i) = TCParams.mintim(1,i);
        TCParams.secondpeaktime(1,i)= TCParams.maxtim(1,i); 
    end
end 

for i = 1:length(TCParams.mintim)
     for p = TCParams.firstpeaktime(1,i):-1:(TCParams.secondpeaktime(1,i)-1) %increments in between peaks
         if((timecourse(p,i) >= 0 &&  timecourse(p-1,i) < 0) || (timecourse(p,i) < 0 &&  timecourse(p-1,i) >= 0)) %if a sign change
             TCParams.zerocrossing(1,i) = p;
             break;
         end
     end
end

TCParams.pk_dot =1 - abs(TCParams.minmaxadd) ./ TCParams.minmaxabsadd; %peak degree of transience 1=biphasic 0 = monophasic 

%put some into time
TCParams.mintim = timevector(TCParams.mintim);
TCParams.maxtim = timevector(TCParams.maxtim);
TCParams.zerocrossing = timevector(TCParams.zerocrossing);
i=1;
  if(TCParams.mintim(1,i) < TCParams.maxtim(1,i)) %if on 
      TCParams.firstpeaktime(1,i) = TCParams.maxtim(1,i);
      TCParams.secondpeaktime(1,i)= TCParams.mintim(1,i); 

  else %if off
      TCParams.firstpeaktime(1,i) = TCParams.mintim(1,i);
      TCParams.secondpeaktime(1,i)= TCParams.maxtim(1,i); 
      
  end
  
end




