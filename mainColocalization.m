%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version of June 21st 2021                                               %
% Authored by Boris Louis,Raffaele Vitale, Susana Rocha and Aline Acke.   %
%                                                                         %
%                                                                         %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
%% User input
%Input single-molecule detection
locROI = 70; %radius in pixel
chi2 = 100; %certainty threshold for initial detection (500)
FWHM = 4; %full width half maximum of the PSF

%Input for simulation to check p-value
nsim=1e4; %number of simulation to be ran

% filename for saving data (excel format)
filename = ' post covid test BRD4 acetylation cell 50 rotation 180.xlsx'; 

%% Data loading
%Allow user to select a folder containing .tif files to analyze
% The folder is expected to contain one of each type of file:
% 1) Nucleus staining (facultative), should contain nucleusStaining in the
% filename
% 2) Particles (mandatory), should contain particle in the filename
% 3) epigenetic marker (mandatory) should contain epigeneticMarker in the
% filename

[folder] = uigetdir();
%get list of files in selected directory
file2Analyze = dir(folder);

%keep only files (remove directories from the list)
file2Analyze([file2Analyze.isdir]) = [];

%keep only compatible files
idx = contains({file2Analyze.name},'.tif','IgnoreCase',true);

file2Analyze = file2Analyze(idx);

%check that we have some files left before we continue
assert(~isempty(file2Analyze),'Could not find compatible files in the given directory, please check the filenames and extensions');

% Data loading Nucleus Staining

nucleusIm = Load.Data(file2Analyze,'nucleusStain');

% Data loading Particles
imParticle = Load.Data(file2Analyze,'particle');
nFrames = size(imParticle,3); 
% check that we found some particle data otherwise throw an error 
assert(~isempty(imParticle),'Could not find particle file with .tif extension, please check the selected directory');
    
% Data loading - DNA marker

imMarker = Load.Data(file2Analyze,'epigeneticMarker');
assert(~isempty(imMarker),'Could not find epigeneticMarker file with .tif extension, please check the selected directory');


%% Segment nucleus (only based on intensity, can be optimized)
% get threshold with help of user
if ~isempty(nucleusIm)
    nFrames = size(nucleusIm,3); 
    fr = round(nFrames/2); %use frame on the middle of the stack
    currentFrame = nucleusIm(:,:,fr);% frame changed to nr 1, change back 1->fr

    [tHold] = expMic.getTh(currentFrame,'Segment Nucleus');
    close(gcf)
    imNuc=uint16(zeros(size(nucleusIm)));% segmented nucleus
    % get number of nucleus
    %user imput: number of cells
    figure, imagesc(currentFrame); axis equal tight
    nr = str2double(inputdlg('How many cells are in the image?', 'Define number of cells', 1));
    close (gcf)

    h = waitbar(0, 'Segmenting nucleus...');
    for i = 1:nFrames
        imNuc(:,:,i) = expMic.segCellNuc(currentFrame,tHold,nr);
        waitbar(i/nFrames, h);
        if i == 1
            imagesc(imNuc(:,:,1))
            pause
        end
    end
    close(h)
    
else
    imNuc = [];
    
end

%% Localization
locPos = cell(nFrames,1);

% is the segmentation is empty then we replace it by a matrix of ones so
% multiplication by it will not do anything.
if isempty(imNuc)
    imNuc = ones(size(imParticle));
end
h = waitbar(0,'Localizing ...');
for i = 1 : nFrames
    currentFrame = double(imParticle(:,:,i).*imNuc(:,:,i));%need double for localization
    %use DAPI to remove particles outside the nucleus, IM_nuc
    
    [pos] = Localization.smDetection(currentFrame,locROI,FWHM,chi2);
    %pos is (y,x) because y is the first dimension in matrices in matlab
    pos = round(pos); %round for ROI
    if ~isempty(pos)
        %Pre- allocation
        SRPos = table(zeros(size(pos, 1), 1), zeros(size(pos, 1), 1), zeros(size(pos, 1), 1), zeros(size(pos, 1), 1),'VariableNames',{'frame','x','y', 'int'});
        
        for j = 1:size(pos,1)
            ROI = currentFrame(pos(j,1)-locROI:pos(j,1)+locROI,pos(j,2)-locROI:pos(j,2)+locROI);

            %Create a grid for fitting
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
figure 
%for i = 1 : nFrames
i=1;%for checking 1 frame only
imagesc(imParticle(:,:,i), [0 15000]);
colormap('jet'); title(['FRAME: ' num2str(i)]);
hold on
try
    plot(locPos{i}.x,locPos{i}.y,'wo');
    %end
    pause
catch e
end

%% eliminate partciles in just one plane + get the brighest plane for each particle

cod = table2array(vertcat(locPos{:})); 

res = Misc.H_autotrack4(cod, 5, 0); %connect detected molecules

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
%selPos = pt;
selPos = array2table(pt, 'VariableNames', {'x', 'y', 'frame', 'particle', 'intensity'});

for i = 1 : length (selPos.frame)
    %fr = selPos(i, 3); %frame where point is the brightest
    fr = selPos.frame(i);
    currentFrame = double(imMarker(:,:,fr)); %frame where particle is brightest
    
    realVal = getIntensities (currentFrame, selPos(i, {'x', 'y'}), 1); %get descriptors from 3x3 area
    realVal = array2table(realVal, 'VariableNames', {'mean', 'median', 'sum', 'max'});
    % simulate nsim point in the same region of interest = nucleus selected
    % by DAPI
    mask = imNuc(:,:, fr); %get DAPI-selected nucleus
    
    
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
    
    %save the mean value of the marker intensity
    selPos.mean_marker(i) = realVal.mean;
    selPos.median_marker(i) = realVal.median;
    

end

%% Get distance to mask

if ~isempty(imNuc)
    % get the centroid position of the cell
    center = regionprops(imNuc,'Centroid');
    center = center(1).Centroid;
    
    %get the countour position
    contour = cell(size(imNuc,3),1);
    for i = 1:size(imNuc,3)
       [pContour] = bwboundaries(imNuc(:,:,i));
       if length(pContour)>1
           error('Found more than one contour please check');
       end
           
       contour{i} = pContour{1};
        
    end
    
    %calculate distance between particle and cell
    for i = 1:height(selPos)
        currContour = contour{selPos.frame(i,:)};
        
        dist2Center = sqrt((selPos.x(i,:)-center(1)).^2 + (selPos.y(i,:)-center(2)).^2);
        
        dist2Membrane = min(sqrt((selPos.x(i,:)-currContour(:,2)).^2 + (selPos.y(i,:)-currContour(:,1)).^2));
        
        selPos.dist2Center(i) = dist2Center;
        selPos.dist2Membrane(i) = dist2Membrane;
                
    end    
end

%% get the Excell file of the p-value
T = selPos; 
writetable(T,filename,'Sheet',1,'Range','A1')