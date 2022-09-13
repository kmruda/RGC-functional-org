
%% Fit spatial RFs to a center/surround difference of Gaussians
% Marija Rudzite

% dataruns: structure with stimulus information and spiking responses of RGCs

% Create spatial tuning curves
for j=1:l_childs %for each light level 
                  
    l_rgb = length(dataruns.(data_names{j}).stimulus.params.RGB);
    l_sp = length(dataruns.(data_names{j}).stimulus.params.SPATIAL_PERIOD);
    SP=dataruns.(data_names{j}).stimulus.params.SPATIAL_PERIOD;       
    SF=(1./SP);                        

    %group spatial period and contrast
    for rgb=1:l_rgb %for each contrast
        for sp=1:l_sp %for each spatial period
            stim_groups_trials{sp,rgb}=[]; 
            for t=1:length(dataruns.(ch_names{j}).stimulus.trials) %for each trial
                if dataruns.(ch_names{j}).stimulus.trials(t).SPATIAL_PERIOD == dataruns.(ch_names{j}).stimulus.params.SPATIAL_PERIOD(sp)...
                        && dataruns.(ch_names{j}).stimulus.trials(t).RGB == dataruns.(ch_names{j}).stimulus.params.RGB{rgb}
                   stim_groups_trials{sp,rgb} = [stim_groups_trials{sp,rgb}(:)',t];
                end
            end
        end
    end

    bs = 0.05;
    samplerate = 1/bs; 

    for i=1:l_type %for each of the types 
        l_cells_in_type=length(IDZ.(ch_names{j}){1,i}); %number of mapped cells within a type 
        get_cell_indices(dataruns.(ch_names{j}), IDZ.(ch_names{j}){1,i});  
     
        All_notnorm1=[];
        All_notnorm2=[];  
        norm_cont1=[];
        norm_cont2=[];
        MAll_contrast_1=[];
        MAll_contrast_2=[];
        All_norm1=[];
        All_norm2=[];
        PS_1=[];
       for k=1:l_cells_in_type %for each mapped cell within a cell type 
            cells = get_cell_indices(dataruns.(ch_names{j}), IDZ.(ch_names{j}){1,i}(1,k)); %counts through the cell indeces  
            for rgb=1:l_rgb %for each contrast
                for sp=1:l_sp %for each spatial period

                    % get PSTH for each cell/contrast/SP
                    PSTH_sp_cells= get_psth(dataruns.(ch_names{j}).spikes{cells},...
                    dataruns.(ch_names{j}).stimulus.triggers(stim_groups_trials{sp,rgb}),...
                    'plot_hist',false,'stop',[8],'bin_size',0.05,'foa',0);

                    [xv,ps,ft]=PowerSpectrumFinder(PSTH_sp_cells,samplerate); %calculate power spectrum 
                    [~,f1_ind]=min(abs(xv-2)); %find xv index closest to 2 Hz
                    ind_array=[f1_ind-1, f1_ind, f1_ind+1];%add next 2 points closest to 2 HZ
                    for b=1:3
                        ps_array(b)=ps(ind_array(b));
                        xv_array(b)=xv(ind_array(b));
                    end 
                    f1area=trapz(xv_array,ps_array); %take integral under the power spectrum @ 2hz peak
                    fall=trapz(xv,ps);
                    PS(sp)=ps(f1_ind);
                    ps_point(sp)=f1area/fall;
                    F1Area(sp)=f1area;
                end
                PS1(rgb,:)=PS; %make a matrix of all data values (per cell) per contrast 
            end            

            PS_1=(PS1(1,:)./max(PS1(1,:))); 
            PS_2=(PS1(2,:)./max(PS1(2,:))); 

            All_norm1(k,:)= PS_1;
            All_norm2(k,:)= PS_2;

            All_notnorm1(k,:)=PS1(1,:); %matrix of all data poiints per cell type with first contrast
            All_notnorm2(k,:)= PS1(2,:); %matrix of all data poiints per cell type with second contrast

            MAll_contrast_1=mean(All_notnorm1,1); %average over cells in same cell type in contrast 1 
            MAll_contrast_2=mean(All_notnorm2,1); %average over cells in same cell type in contrast 2 

            norm_cont1=MAll_contrast_1./max(MAll_contrast_1); %average data normalization contrast 1
            norm_cont2=MAll_contrast_2./max(MAll_contrast_2); %average data normalization contrast 2

            fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{1}).alldata= All_notnorm1; %all data contrast 1, not normalized
            fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{2}).alldata= All_notnorm2; %all data contrast 2, not normalized 
            fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{1}).normalldata= All_norm1; %all data contrast 1, normalized 
            fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{2}).normalldata= All_norm2; %all data contrast 2, normalized 

        end
       fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{1}).data=norm_cont1; %averaged, normalized data for contrast 1  
       fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{2}).data=norm_cont2; %averaged, normalized data for contrast 2 
       fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{1}).avenotnorm = MAll_contrast_1;%averaged, not normalized data for contrast 1  
       fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{2}).avenotnorm=MAll_contrast_2;%averaged,not normalized data for contrast 2 
    end
 
end

%Make a cosine matrix 
pi_div = 24;
phase_range = pi/pi_div:pi/pi_div:2*pi;
num_phase = length(phase_range);
for sf=1:length(SF)             
    for phase_index=1:num_phase
        x=0:1:699;
        y_term= cos(SF(sf)*x+ phase_range(phase_index));
        y_mat = repmat(y_term, 700, 1);
        cos_term(:,phase_index)=y_mat(:);
    end        
    All_cos_term(:,:,sf)=cos_term;
end

% Fit tuning curves to difference of Gaussians
for rgb=1:l_rgb %for each contrast    
    for j=1:l_childs % for each light level 
        for i=1:l_type %for each cell type
            l=size(fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{rgb}).alldata);
            for k=1:l(1) %for each cell 
            
                SP=dataruns.(data_names{j}).stimulus.params.SPATIAL_PERIOD;
                SF=(1./SP);
                x_data=SF.*2*pi;           
                y_data=fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{rgb}).normalldata(k,:);
                b0 = [45, 100, 5, 1]; %setting initial parameters
                LB= [5,10,0.000001,0];
                UB= [150,350,10,10];         

                %FMINCON
                %contrain equations
                A1 = [1.5, -1, 0, 0];
                A2 = [-10, 1, 0, 0];

                myfun = @(b)get_DOG_error(b,x_data,y_data,All_cos_term); 
                problem = createOptimProblem('fmincon','objective',myfun,'x0',b0,'Aineq',[A1;A2],'bineq',[0,0],'lb',LB,'ub',UB);         
                gs=GlobalSearch('Display','iter','NumTrialPoints',250,'StartPointsToRun','bounds');
                rng(14,'twister'); %global search going through different parameters until smallest difference is found  
                [AMRb,fval]=run(gs,problem);
                fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{rgb}).fit_ALL(k,:)=AMRb;
                fmin2Results.(ch_names{j}).(types_interest{i}).(RGB{rgb}).lsqer(k,:)=fval;

            end
        end
    end
end

%%
function DOG_error = get_DOG_error(b,SF,y_data,All_cos_term)
    %Fitting function
    temp_params.center_radius = b(1);
    temp_params.surround_radius = b(2);
    temp_params.center_scale = b(3);
    temp_params.surround_scale = b(4);
    temp_params.x_size = 700;
    temp_params.y_size = 700;
    temp_params.center=[350,350]; 
    temp_rf = make_gaussian(temp_params);  

    num_phase = size(All_cos_term,2);
    DOG_term=repmat(temp_rf(:),1,num_phase);
    ALL_DOG_term=repmat(DOG_term,1,1,length(SF));

    dot_product_AMR=dot(All_cos_term,ALL_DOG_term);
    tuning_AMR=max(dot_product_AMR);
    sf_tuning=reshape(tuning_AMR,[length(SF),1]);
    sf_tuning = (sf_tuning')./max(sf_tuning);

    DOG_error=sqrt(mean((sf_tuning-y_data).^2));
end
     
