function [bblist] = global_detection(I,submean,meanI,extractPara,owndetector)

%==========================================================================
% Global detection -> detect the whole single frame 
%==========================================================================
% Demo version v1.1
% 2022-10-16
% © Xicheng Wang, Yun Feng, Dmitry Grishchenko, Pavel Kudinov
%==========================================================================

%==========================================================================
% 1. Bubble detection
%==========================================================================
% Demo version v1.0
% 2022-06-17
% © Yun Feng, Xicheng Wang, Dmitry Grishchenko, Pavel Kudinov
%==========================================================================
    
% Subtract mean of all images to remove the background of each image
%It improves results of SBTV
if submean == 1 
    I = mat2gray(double(I)-meanI);  % subtract mean value to remove noisy background
else
    I = mat2gray(I);  % subtract mean value to remove noisy background
end

% Individual detection using  (grid offset + adaptive threshold)
% |area| |center X| |center Y| |Bbox: x0 y0 L W|
[bblist,I] =  griddetection(I,extractPara); 
end

