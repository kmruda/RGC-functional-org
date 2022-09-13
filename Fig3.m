
%% Perform SVD on 3-D RFs

% avg_stas: structure of averaged STAs for each cell type/light level

types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
ndf_names = {'NDF3';'NDF0'};

height = 40; width = 53;
for type=1:length(types_interest)
    for j=1:length(ndf_names)
        this_sta = avg_stas.(ndf_names{j}){type};
        sta3d = squeeze(this_sta(:,:,1,:));
        sta2d = reshape(sta3d,height*width,[])';
        [usta,ssta,vsta] = svd(sta2d);
        temp = diag(ssta);
        svspec.(types_interest{type})(j,:) = temp./sum(temp); % explained variance of full STA
    end
end

% plot variance explained by rank 1 RF
figure;
for j = 1:length(ndf_names)
    for type = 1:length(types_interest)
        xs = j+rand(1)./3;
        plot(xs,svspec.(types_interest{type})(j,1),'o');hold on
    end
end
xlim([0.5 2.5])
