function mask = segCellNuc(IM,tHold,nr)
%SEGCELLNUC segment the cell nucleus
%   Detailed explanation goes here
%   nr is number of nucleus in the image

% apply treshold
BWpre = IM>tHold;
% clear border
BWpre = imclearborder(BWpre);
% remove small things
BWpre = bwareaopen(BWpre,5);
% smooth
se = strel('disk',5);
BWpre = imclose(BWpre,se);
BWpre = imopen(BWpre,se);
% remove small things
BWpre = bwareaopen(BWpre,100);
BWpre = ~BWpre;
BWpre = bwareaopen(BWpre,100);
BWpre = ~BWpre;

% get boundaries
[B,~] = bwboundaries(BWpre);
% size of the boundries
Blength = cellfun(@length,B);

% cells must be the largest of them
[~, idx] = maxk(Blength,nr);
B = B(idx);

%smooth
for i = 1:nr
    B{i} = SDcalc.smoothBoundary( B{i} );
end

Ball = round(cell2mat(B));
mask = accumarray([Ball(:,1),Ball(:,2)],1,size(IM));
mask = imfill(mask);

end
