clear
clc
close all
%% User input
locROI = 6; %radius in pixel
chi2 = 80; %certainty threshold for initial detection
FWHM = 2; %full width half maximum of the PSF
%% Data loading
%open UI to search file
[fileName,folder,~] = uigetfile('*.tif','Pick a tif file to analyze');
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
locPos = cell(nFrames);
h = waitbar(0,'Localizing ...');
for i = 1 : nFrames
    currentFrame = double(IM(:,:,1));%need double for localization
    
    [pos] = Localization.smDetection(currentFrame,locROI,FWHM,chi2);
    %pos is (y,x) because y is the first dimension in matrices in matlab
    pos = round(pos); %round for ROI
    SRPos = table(pos(:,2),pos(:,1),'VariableNames',{'x','y'});
    
    for j = 1:size(pos,1)
        ROI = currentFrame(pos(j,1)-locROI:pos(j,1)+locROI,pos(j,2)-locROI:pos(j,2)+locROI);
        
        [y,x,e] = Localization.phasor(ROI);
        
        if e>1.5 || e<0.6
            warning('Molecules seems very elongated, something might be wrong there,perhaps out of focus');
        end
        
        %recalculate xy position in original coordinate system;
        SRPos.x(j) = x+pos(j,2);
        SRPos.y(j) = y+pos(j,1);
        
    end
    
    locPos{i} = SRPos;
    
    waitbar(i/nFrames,h,sprintf('Localizing Frame: %d / %d',i, nFrames));
end
close(h)
%% plot to check
figure 
imagesc(IM(:,:,1));
colormap('hot');
hold on
plot(locPos{1}.x,locPos{1}.y,'wo');
