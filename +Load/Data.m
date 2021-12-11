function [data] = Data(fileList,keyword)

    idx = contains({fileList.name},keyword,'IgnoreCase',true);

    switch(sum(idx))
    case 0

        disp(['No ' keyword ' file found']);
        nucleusIm = [];
    
    case 1

        fullPath = [fileList(idx).folder filesep fileList(idx).name]; 
        
    otherwise
        error('More than one nucleus staining found, this cannot be handled at the moment');
    end
    
    %get the info of the file
    [fileInfo] = Load.Movie.tif.getinfo(fullPath);

    %get the total number of frame in the file and load them
    frames = 1:fileInfo.Frame_n;
    nFrames = length(frames);
    data = Load.Movie.tif.getframes(fullPath,frames);

end