clear
clc
close all

%% Data loading DAPI
%open UI to search file
[fileName,folder,~] = uigetfile('*DAPI*.tif','Pick a tif file to analyze');
cd(folder)
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);

%check the extension is correct
assert(strcmp(ext,'.tif'),'The file picked is not a tif file');

%get the info of the file
[fileInfo] = Load.Movie.tif.getinfo(fullPath);

%get the total number of frame in the file
frames = 1:fileInfo.Frame_n;
nFrames = length(frames);
IM_DAPI = Load.Movie.tif.getframes(fullPath,frames);

%% Segment nucleus (only based on intensity, can be optimized)
% get threshold with help of user
fr = round(nFrames/2); %use frame on the middle of the stack
currentFrame = IM_DAPI(:,:,fr);

[tHold] = expMic.getTh(currentFrame,'Segment Nucleus');
close(gcf)
IM_nuc=uint16(zeros(size(IM_DAPI)));% segmented nucleus
% get number of nucleus
%user imput: number of cells
figure, imagesc(currentFrame); axis equal tight
nr = str2double(inputdlg('How many cells are in the image?', 'Define number of cells', 1));
close (gcf)

h = waitbar(0, 'Segmenting nucleus...');
for i = 1:nFrames
    IM_nuc(:,:,i) = expMic.segCellNuc(currentFrame,tHold,nr);
    waitbar(i/nFrames, h);
end
close(h)


%% User input
locROI = 20; %radius in pixel
chi2 = 500; %certainty threshold for initial detection
FWHM = 4; %full width half maximum of the PSF
%% Data loading VIRUS
%open UI to search file
[fileName,folder,~] = uigetfile('*irus*.tif','Pick a tif file to analyze');
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);

%check the extension is correct
assert(strcmp(ext,'.tif'),'The file picked is not a tif file');

%get the info of the file
[fileInfo] = Load.Movie.tif.getinfo(fullPath);

%get the total number of frame in the file
frames = 1:fileInfo.Frame_n;
nFrames = length(frames);
IM = Load.Movie.tif.getframes(fullPath,frames);

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
%% plot to check
% figure 
% for i = 1 : nFrames
% imagesc(IM(:,:,i), [0 15000]);
% colormap('jet'); title(['FRAME: ' num2str(i)]);
% hold on
% try
% plot(locPos{i}.x,locPos{i}.y,'wo');
% end
% pause
% end

%% eliminate partciles in just one plane + get the brighest plane for each particle

cod = table2array(vertcat(locPos{:})); 

res = H_autotrack4(cod, 3, 0); %connect detected molecules

pt = []; %points to be analysed

for i = 1:max(res(:,4))
    tr = res(res(:,4)==i, :);
    if size(tr, 1) > 1 % particle is detected in 2 consecutive planes
        [ind, ~]=find(tr==max(tr(:,5)));
        pt = cat(1, pt, tr(ind, :));
    end
end

%pt = array2table(pt, 'VariableNames', {'x', 'y', 'frame', 'particleNr', 'int'});


%% Data loading - DNA marker
%open UI to search file
[fileName,folder,~] = uigetfile('*arker*.tif','Pick a tif file to analyze');
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);

%check the extension is correct
assert(strcmp(ext,'.tif'),'The file picked is not a tif file');

%get the info of the file
[fileInfo] = Load.Movie.tif.getinfo(fullPath);

%get the total number of frame in the file
frames = 1:fileInfo.Frame_n;
nFrames = length(frames);
IM_marker = Load.Movie.tif.getframes(fullPath,frames);

%% calculate p-value of selected points
%selPos = pt;
selPos = array2table(pt, 'VariableNames', {'x', 'y', 'frame', 'particle', 'intensity'});

%pval = cell2table(cell(0,4), 'VariableNames', {'mean', 'max', 'sum', 'median'});

for i = 1 : length (selPos.frame);
    %fr = selPos(i, 3); %frame where point is the brightest
    fr = selPos.frame(i);
    currentFrame = double(IM_marker(:,:,fr)); %frame where particle is brightest
    %currentFrame = double(IM(:,:,fr)); %load virus image for control
    %currentFrame = double(max(IM(:,:,fr), IM_marker(:,:,fr))); %load virus image for control + marker
    
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

%% get the Excell file of the p-value

T = selPos;
filename = 'MLV WT cell 1 sample 4.4.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1')
