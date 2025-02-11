function [bbinfo] = localRoiDdetection(image,ROI,extractPara)
% localRoiDdetection
% Detect bubbles from local ROI with classical detector and Faster RCNN detector
%##########################################################################
%--------------------------------------------------------------------------
% 'image'        -> grey scale figure
% 'extractPara'  -> object stores user defined parameters for extraction
%--------------------------------------------------------------------------

% Avoid blackground
if ~isnan(extractPara.trsh)
    
    % Bubble detection with classical detector
    rp=regionprops('table', image<extractPara.trsh); 
    % Filter the bubble outside the area range
    szfilt = ([rp.Area] < extractPara.areafilt(1))|([rp.Area] > extractPara.areafilt(2)); 
    rp(szfilt,:) = [];  
    bbinfo = table2array(rp);
    
    if ~isempty(bbinfo)
        % add mask
        for i=1:length(bbinfo(:,1))
            image = addMask(image,bbinfo(i,[2:3,6:7]),extractPara.maskMargin,1);
        end
        % move bubble center to global coordinate
        bbinfo(:,[2:3]) = [ROI(1),ROI(2)] + rp.Centroid - [1 1];
        bbinfo(:,[4:5]) =[ROI(1),ROI(2)]+rp.BoundingBox(:,1:2)-[1 1];
    end

else
    bbinfo = [];
end

end
