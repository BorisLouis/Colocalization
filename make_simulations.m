clear
clc
close all


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

%% calculate similate data with variable amount of overlap
%selPos = pt;
toadd = inputdlg('How much signal do you want to add?', 'Add virus signal to marker channel', 1, {'0.1'});
outfile = [fullPath(1:findstr(fullPath, '.lif')-1) '_sim' toadd{1} '.tif'];
toadd = str2double(toadd{1}); %fraction of virus signal to be added

%for the first frame
currentFrame = uint16(double(IM_marker(:,:,1)) + toadd.*double(IM(:,:,1)));
imwrite(currentFrame, outfile)
h = waitbar(0, 'Calculating data...');
 
for i = 2 : nFrames
    waitbar(i/nFrames, h);
   currentFrame = uint16(double(IM_marker(:,:,i)) + toadd.*double(IM(:,:,i))); %frame where particle is brightest
   imwrite(currentFrame, outfile, 'WriteMode','append');
end
close(h)