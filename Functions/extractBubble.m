function [rp] = extractBubble(image,ROI,extractPara)
%extractBubble Summary of this function goes here
%   Detect bubbles from image based on local adpative threshold 
%##########################################################################
%--------------------------------------------------------------------------
% 'image'        -> grey scale figure
% 'extractPara'  -> object stores user defined parameters for extraction
%--------------------------------------------------------------------------
% Avoid blackground
if ~isnan(extractPara.trsh)
    % Bubble extraction 
    rp=regionprops('table', image<extractPara.trsh); 
    % Filter the bubble outside the area range
    szfilt = ([rp.Area] < extractPara.areafilt(1))|([rp.Area] > extractPara.areafilt(2)); 
    rp(szfilt,:) = [];  
    
    if ~isempty(rp) 
        % move bubble center to global coordinate
        rp.Centroid = [ROI(1),ROI(2)] + rp.Centroid - [1 1];
        rp.BoundingBox(:,1:2) =[ROI(1),ROI(2)]+rp.BoundingBox(:,1:2)-[1 1];
    end

else
    rp = [];
end
end
