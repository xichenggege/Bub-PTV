function [ROI] = createROIwindow(pos,height,width,im_H,im_W)
%CREATEROIWINDOW Summary of this function goes here
%   extend ROI window based on width and length
%--------------------------------------------------------------------------
% left top corner and right bottom corner
% make sure ROI is in the image
X0 = round(pos(1) - width);
Y0 = round(pos(2) - height);
X1 = X0 + 2*width - 1;
Y1 = Y0 + 2*height - 1;
% not exceed to figure
if Y0<1
    Y0 = 1;
end

if Y1>im_H
    Y1 = im_H;
end

if X1>im_W
    X1 = im_W;
end

 if X0<1
    X0 = 1;
end
ROI = [X0 Y0 X1 Y1];

end

