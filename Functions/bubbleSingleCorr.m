function [corr] = bubbleSingleCorr(bb1,bb2,error,corrPara)
%bubbleSingleCorr correlate a pair of bubbles of two images
% The bubble in left would be registered first and then correlated by
% epipolar constraint

% No transformation as fMatrix from calibration with unregistered images

% Corresponding epipolar line in left 
epiLines = epipolarLine(corrPara.fMatrix',bb2);
% If it satify epipolar constraint
h = (bb1(:,1)*epiLines(1)+epiLines(3))/(-epiLines(2)); % height in image 1
if abs(h-bb1(:,2))<= error
    corr = 1;
else
    corr = 0;
end
corr = logical(corr);
end




