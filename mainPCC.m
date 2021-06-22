clear
clc

%%
[fileName,folder,~] = uigetfile('*DAPI*.tif','Pick a tif file to analyze');
cd(folder)
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);
fileDAPI = fileName;

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
    if i == 1
        imagesc(IM_nuc(:,:,1))
        pause
    end
end
close(h)

%% Data loading VIRUS
%open UI to search file
[fileName,folder,~] = uigetfile('*irus*.tif','Pick a tif file to analyze');
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);
fileVirus = fileName;

%check the extension is correct
assert(strcmp(ext,'.tif'),'The file picked is not a tif file');

%get the info of the file
[fileInfo] = Load.Movie.tif.getinfo(fullPath);

%get the total number of frame in the file
frames = 1:fileInfo.Frame_n;
nFrames = length(frames);
IM = Load.Movie.tif.getframes(fullPath,frames);

%% Data loading - DNA marker
%open UI to search file
[fileName,folder,~] = uigetfile('*arker*.tif','Pick a tif file to analyze');
fullPath = [folder filesep fileName]; 
[~,~,ext] = fileparts(fullPath);
fileMarker = fileName;

%check the extension is correct
assert(strcmp(ext,'.tif'),'The file picked is not a tif file');

%get the info of the file
[fileInfo] = Load.Movie.tif.getinfo(fullPath);

%get the total number of frame in the file
frames = 1:fileInfo.Frame_n;
nFrames = length(frames);
IM_marker = Load.Movie.tif.getframes(fullPath,frames);


%% CALCULATE PCC
PCC = zeros(size(IM,3),1);
for fr = 1:size(IM,3)
    sel_IM = IM(:,:,fr); sel_IM = sel_IM(:); sel_IM (IM_nuc(:,:,fr) == 0) = []; %delete regions outide nucleus 
    sel_IM_marker = IM_marker(:,:,fr); sel_IM_marker = sel_IM_marker(:); sel_IM_marker (IM_nuc(:,:,fr) == 0) = []; %delete regions outide nucleus 
    PCC(fr) = corr2(sel_IM, sel_IM_marker);
end

PCC
mean(PCC)
 
%%
save(['Analysis PCC ' fileVirus '.mat'], 'PCC', 'fileMarker', 'fileDAPI', 'fileVirus', 'IM_nuc', 'IM_DAPI', 'IM_marker', 'IM');


