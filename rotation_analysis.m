clear
clc
close all

%% Data loading 
%open UI to search file
[fileName,folder,~] = uigetfile('Analysis*.mat','Pick a Analysis file');
cd(folder)
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);
file = fileName;
load(file)

IMori=IM; % save unrotated/original virus image
nFrames=size(IM, 3);

%% User input for point detection
locROI = 20; %radius in pixel
chi2 = 70; %certainty threshold for initial detection
FWHM = 4; %full width half maximum of the PSF

selPos_all = cell(1,4);
%% Calculate position and p-values for 4 different rotations
for r=0:90:270 %check angles 0, 90, 180, 270
    IM=imrotate(IMori, r, 'nearest', 'crop');
    %% Localization
    locPos = cell(nFrames,1);
    h = waitbar(0,'Localizing ...');
    for i = 1 : nFrames
        currentFrame = double(IM(:,:,i).*IM_nuc(:,:,i));%need double for localization
        %use DAPI to remove particles outside the nucleus, IM_nuc
        
        [pos] = Localization.smDetection(currentFrame,locROI,FWHM,chi2);
        %pos is (y,x) because y is the first dimension in matrices in matlab
        pos = round(pos); %round for ROI
        if ~isempty(pos)
            SRPos = table(zeros(size(pos, 1), 1), zeros(size(pos, 1), 1), zeros(size(pos, 1), 1), zeros(size(pos, 1), 1),'VariableNames',{'frame','x','y', 'int'});
            
            for j = 1:size(pos,1)
                ROI = currentFrame(pos(j,1)-locROI:pos(j,1)+locROI,pos(j,2)-locROI:pos(j,2)+locROI);
                
                
                [X,Y] = meshgrid(pos(j,2)-locROI:pos(j,2)+locROI,pos(j,1)-locROI:pos(j,1)+locROI);
                domain(:,:,1) = X;
                domain(:,:,2) = Y;
                
                [gPar] = Localization.Gauss.MultipleFitting(ROI,pos(j,2),pos(j,1),domain,1);%data,x0,y0,domain,nbOfFit
                
                SRPos.frame(j) = i;
                SRPos.x(j) = gPar(5);
                SRPos.y(j) = gPar(6);
                SRPos.int(j) = gPar(1);
                
                
            end
            
            locPos{i} = SRPos;
        end
        
        waitbar(i/nFrames,h,sprintf('Localizing Frame: %d / %d',i, nFrames));
    end
    close(h)
    %% eliminate partiles in just one plane + get the brighest plane for each particle
    
    cod = table2array(vertcat(locPos{:}));
    
    res = H_autotrack4(cod, 3, 0); %connect detected molecules
    
    pt = []; %points to be analysed
    
   for i = 1:max(res(:,4))
    tr = res(res(:,4)==i, :);
    %if size(tr, 1) > 1 % particle is detected in 2 consecutive planes
        [ind, ~]=find(tr==max(tr(:,5)));
        pt = cat(1, pt, tr(ind, :));
    %end
end
%pt = [cod(:,2:3) cod(:,1) ones(size(cod,1), 1) cod(:,4)];
pt = sortrows(pt, 5);
pt = pt(floor(size(pt, 1)*0.1)+1:end, :);

%pt = array2table(pt, 'VariableNames', {'x', 'y', 'frame', 'particleNr', 'int'});
    %% calculate p-value of selected points
    selPos = array2table(pt, 'VariableNames', {'x', 'y', 'frame', 'particle', 'intensity'});
    
    for i = 1 : length (selPos.frame)
        %fr = selPos(i, 3); %frame where point is the brightest
        fr = selPos.frame(i);
        currentFrame = double(IM_marker(:,:,fr)); %frame where particle is brightest
        
        realVal = getIntensities (currentFrame, selPos(i, {'x', 'y'}), 1); %get descriptors from 3x3 area
        realVal = array2table(realVal, 'VariableNames', {'mean', 'median', 'sum', 'max'});
        % simulate nsim point in the same region of interest = nucleus selected
        % by DAPI
        mask = IM_nuc(:,:, fr); %get DAPI-selected nucleus
        
        nsim=1e4;
        index  = find(mask);
        select = index(randperm(length(index), nsim));
        [rand(:,2), rand(:,1)] = ind2sub(size(mask), select); %10000 random points inside the DAPI-defined nucleus
        
        simVal = arrayfun(@(x) getIntensities(currentFrame, rand(x,:), 1), 1:size(rand, 1), 'UniformOutput', false);
        simVal = reshape(cell2mat(simVal), 4, [])';
        simVal = array2table(simVal, 'VariableNames', {'mean', 'median', 'sum', 'max'});
        
        %calculate pValue for point i
        selPos.pval_mean(i) = (sum(realVal.mean<simVal.mean) + 1)/(nsim+1);
        selPos.pval_median(i) = (sum(realVal.median<simVal.median) + 1)/(nsim+1);
        selPos.pval_max(i) = (sum(realVal.max<simVal.max) + 1)/(nsim+1);
        selPos.pval_sum(i) = (sum(realVal.sum<simVal.sum) + 1)/(nsim+1);
        
    end
    selPos_all{r/90+1}=selPos;
    disp('rotating image...');
end

%% calculate % p-values for diffenrent conditions
res=[];
sz = reshape(cell2mat(cellfun(@(x) size(x), selPos_all, 'UniformOutput', false)), [2,4]);

res(1,:) = sz(1,:);
res(2,:) = cell2mat(arrayfun(@(x) length(find(selPos_all{1,x}.pval_mean<0.05)), [1:4], 'UniformOutput',false));
res(3,:)= res(2,:)./res(1,:).*100;
pval = array2table(res', 'VariableNames', {'detec_pt', 'pt<0.05', '%'});
disp(pval)
save(['Rotation ' fileVirus '.mat'], 'res', 'pval', 'selPos_all', 'fileMarker', 'fileDAPI', 'fileVirus', 'IM_nuc', 'IM_DAPI', 'IM_marker', 'IM');