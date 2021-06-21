function [intVal] = getIntensities (im, cod, n)
%function that get the mean, median, sum, max intensity of a area n
%around point in cod (x,y), in image im

% if n=1 ,it is an area of 3x3
% if n=2, it is an area of 5x5

%structure of intVal
% 1 = mean; 2 = median; 3 = sum; 4 = max;

if istable(cod) || isstruct(cod)
    cod2=[cod.x cod.y];
    cod=round(cod2);  
else
    cod=round(cod(:,1:2));
end

try
    ROI = im(cod(2)-n:cod(2)+n,cod(1)-n:cod(1)+n);
    ROI = ROI(:);
    
    intVal = [mean(ROI) median(ROI) sum(ROI) max(ROI)];
catch
    intVal = [-1 -1 -1 -1];
end


