
%% Quantify extent of space-time tradeoff 

% avg_tt0_per_ret: structure of time for temporal RF to cross the x-axis for each cell type/light level
% avg_radii_per_ret: structure of spatial RF sizes (center radius) for each cell type/light level
% dr_paths: data per retina

types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
all_ndf_names = {'NDF5';'NDF4';'NDF3';'NDF2';'NDF1';'NDF0'};

% normalize x and y values based on average of all types
normerradii = [];normertt0 = []; 
for dr = 1:length(dr_paths)
    for j = 1:length(all_ndf_names)
        if dr <= length(avg_tt0_per_ret.(types_interest{type}).(all_ndf_names{j}))
            fornormx = []; fornormy = [];
            for type = 1:length(types_interest)
                fornormx(type) = avg_radii_per_ret.(types_interest{type}).(all_ndf_names{j})(dr);
                fornormy(type) = -1*avg_tt0_per_ret.(types_interest{type}).(all_ndf_names{j})(dr);
            end
            
            normerx = mean(fornormx);
            normery = mean(fornormy);

            for type = 1:length(types_interest)
                x = avg_radii_per_ret.(types_interest{type}).(all_ndf_names{j})(dr)/normerx;
                y = -1*avg_tt0_per_ret.(types_interest{type}).(all_ndf_names{j})(dr)/normery;

                normerradii{dr}(type,j) = x;
                normertt0{dr}(type,j) = y;            
            end
        end
    end
end
for dr = 1:length(dr_paths)
    for j = 1:length(all_ndf_names)
        if dr <= length(avg_tt0_per_ret.(types_interest{type}).(all_ndf_names{j}))
            for type = 1:length(types_interest)
                normerradii{dr}(normerradii{dr}==0) = nan;
                normertt0{dr}(normertt0{dr}==0) = nan;
            end
        end
    end
end

% fit line for all the data (constrained to 1 slope)
figure;
for dr = 1:length(dr_paths)
    subplot(2,3,dr);
    for j = 1:length(all_ndf_names)
        for type = 1:length(types_interest)
            x = normerradii{dr}(type,j);
            y = normertt0{dr}(type,j);
            plot(x,y,'o','Color',cs_ll(j,:));hold on           
        end
    end
    xlabel('rf radius');ylabel('zero cross');
    
    tempx = normerradii{dr}(:);
    tempy = normertt0{dr}(:);
    [coeffs,S] = polyfit(tempx(~isnan(tempx)),tempy(~isnan(tempy)),1);
    fitX = linspace(0.5,2,100);
    fitY = polyval(coeffs,fitX,S);
    plot(fitX,fitY,'--k')

    temp = corrcoef(tempx(~isnan(tempx)),tempy(~isnan(tempy)));
    r2 = temp(1,2)^2;
    one_slope_r2(dr) = r2;
    
    % fit line for each light level separately
    for j = 1:length(all_ndf_names)
        if ~isnan(normerradii{dr}(1,j))
            tempx = normerradii{dr}(:,j);
            tempy = normertt0{dr}(:,j);
            [coeffs,S] = polyfit(tempx(~isnan(tempx)),tempy(~isnan(tempy)),1);
            fitX = linspace(0.5,2,100);
            fitY = polyval(coeffs,fitX,S);
            plot(fitX,fitY,':','Color',[0.5 0.5 0.5])

            temp = corrcoef(tempx(~isnan(tempx)),tempy(~isnan(tempy)));
            r2 = temp(1,2)^2;
            ll_slope_r2(dr,j) = r2;
        else
            ll_slope_r2(dr,j) = nan;
        end
    end
end

% restrict to light levels where the flexible slope is trusted (better than the fixed slope)
figure;k = 1;
bar(1:length(dr_paths),one_slope_r2,0.25,'FaceColor',[0 0 0]);hold on
for dr = 1:length(dr_paths)
    temp = ll_slope_r2(dr,:);
    temp(temp<one_slope_r2(dr)) = nan; 
    ll_avg = nanmean(temp);
    bar(k+0.4,ll_avg,0.25,'FaceColor',[0.5 0.5 0.5]); 
    k  = k+1;
end
ylabel('r2');

