function [ frameCell, movieInfo, tfl ] = getInfo( path2file )
%GETINFO receives as input the path to the file and gives back all
%infromation about the frames, the movie and total number of frames
%   Detailed explanation goes here
warning('off','all')
tObj = Tiff(path2file,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
idx = strfind(path2file,filesep);
movieInfo.Path   = path2file(1:(idx(end))-1);

assert(tObj.currentDirectory == 1)
header = tObj.getTag(270);

% test here for having multiple image files
imFileStart = strfind(header, '<Image');
imFileEnd   = strfind(header, '</Image>');
nImFiles = length(imFileStart);
assert(nImFiles == length(imFileEnd), 'Problems with image files, size does not match')

if nImFiles > 1
    movieInfo.isMultiImage = true;
else
    movieInfo.isMultiImage = false;
end

frameCell = cell(nImFiles,1);

for imIdx = 1:nImFiles

    % get frame header
    frameHead = header(imFileStart(imIdx):imFileEnd(imIdx)+7);
    % Index frame header to find information relative to each tif index
    [k1, k2, k3, k4, nFrames] = indexFrameHeader(frameHead);
    % init frame info structure
    frameInfo = initFrameInfoStruc(nFrames);
    % fill in frame info
    for i = 1:nFrames
        str1 = frameHead(k1(i):k2(i)-1);
        str2 = frameHead(k3(i):k4(i));
        [ frameInfo(i).C, frameInfo(i).T, frameInfo(i).Z, frameInfo(i).IFD,...
            frameInfo(i).P, frameInfo(i).File, frameInfo(i).Pos,...
            frameInfo(i).expT ] = getInfoFromString( str1, str2 );
    end
    
    frameCell{imIdx} = frameInfo;
end


%Add extrainfo to the movie, in particular, info about camera, max Frame,
%and zStack into movieInfo
movieInfo.isZStack = size(unique({frameInfo.T}),2)==1;
movieInfo.Cam      = str2double(unique({frameInfo.C}));
movieInfo.expT = frameInfo(1).expT;
switch movieInfo.isZStack
    case 0
        for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).T}))+1;
            movieInfo.maxFrame(i) = maxFrame;
        end
    case 1
         for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).Z}))+1;
            movieInfo.maxFrame(i) = maxFrame;
        end
end


tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end
    tObj.nextDirectory();
end

tObj.setDirectory(1)
warning('on','all')

tObj.close

% this is so we dont break backwards compability, but we might want to
% change this into something cleaner later
if ~movieInfo.isMultiImage
    % then we should behave as before, thus we overwrite the cell by the
    % structure.
    frameCell = frameInfo;
    
end

end



function [k1, k2, k3, k4, nFrames] = indexFrameHeader(frameHeader)
% helper function that indexes over the frame header
    k1 = strfind(frameHeader, '<TiffData');
    k2 = strfind(frameHeader, '</TiffData>');
    k3 = strfind(frameHeader, '<Plane');
    k4 = strfind(frameHeader, '"/>');
    k4(k4<min(k3)) = [];
    nFrames = size(k1,2);
    assert(all([length(k2), length(k3), length(k4)]==nFrames),'Plane and tif information do not match');

end

function out = initFrameInfoStruc(nFrames)
% helper function to init an empty frameInfo structure
    out(nFrames).C = [];
    out(nFrames).T = [];
    out(nFrames).Z = [];
    out(nFrames).IFD = [];
    out(nFrames).P = [];
    out(nFrames).File = [];
    out(nFrames).Pos = [];
    out(nFrames).expT = [];
end
