function [bbinfo,I] = griddetection(I,extractPara)
%griddetection Summary of this function goes here
% Segement whole picture and then detect bubbles 
% based on local adpative threshold 
% horizontal and vertical offset to detect bubbles on the grid boundary
%##########################################################################
%--------------------------------------------------------------------------
% 'image'        -> gray scale figure
% 'extractPara'  -> object stores user defined parameters for extraction
%--------------------------------------------------------------------------
% number of segment in length
[im_W,im_L] = size(I);
numSegL     = ceil(im_L/extractPara.segL); % number of segementation in length
numSegW     = ceil(im_W/extractPara.segW); % number of segementation in width (height)
offsetL     = ceil(0.5*extractPara.segL);  % offset in length
offsetW     = ceil(0.5*extractPara.segW);
dynTrsh     = zeros(numSegL,numSegW);  % save dynamic threshold 

%% Global detection to find as much as possible bubbles
bbinfo = [];
% tailor image with ROI
ROI = [1,1,im_L,im_W];
% local adapative thershold value based on 'cannoy' edge detector
extractPara.trsh = adapCannyTrsh(I,0.85,0.95);
% global detection and remove detected bubbles
rp     = extractBubble(I,ROI,extractPara);      
bbinfo = table2array(rp);
% add mask
for i=1:length(bbinfo(:,1))
    I = addMask(I,bbinfo(i,[2:3,6:7]),extractPara.maskMargin,1);
end

%% grid detection, first cycle
bbinfo_cyc1 = [];    
for i=1:numSegL       
    for j=1:numSegW
        % Segemente the whole picture        
        ROI(1) = 1 + (i-1)*extractPara.segL + extractPara.thick;
        ROI(3) = ROI(1) + extractPara.segL - 2*extractPara.thick;
        if ROI(3)>=im_L
            ROI(3) = im_L;
        end

        ROI(2) = 1 + (j-1)*extractPara.segW + extractPara.thick;
        ROI(4) = ROI(2)+extractPara.segW - 2*extractPara.thick;
        if ROI(4)>=im_W
            ROI(4) = im_W;
        end

        % tailor image with ROI
        im_tailor        = I(ROI(2):ROI(4),ROI(1):ROI(3));     
        extractPara.trsh = adapCannyTrsh(im_tailor,0.85,0.95);
        dynTrsh(i,j)     = extractPara.trsh; % save for later offset grid
        rp = extractBubble(im_tailor,ROI,extractPara);
        % collect bubble infomation
        if ~isempty(rp) 
            % area [1], x y [2~3], bounding box [4~7];
            bbinfo_cyc1 = [bbinfo_cyc1;table2array(rp)]; % add to total list
        end
    end
end
% add mask
if ~isempty(bbinfo_cyc1)
    for i=1:length(bbinfo_cyc1(:,1))
        I  = addMask(I,bbinfo_cyc1(i,[2:3,6:7]),extractPara.maskMargin,1); 
    end
end

%% grid detection, second cycle, horizontal offset
bbinfo_cyc2 = [];    
for i=1:numSegL       
    for j=1:numSegW
        % Segemente the whole picture        
        ROI(1) = 1 + (i-1)*extractPara.segL + extractPara.thick + offsetL;
        ROI(3) = ROI(1) + extractPara.segL - 2*extractPara.thick + offsetL;
        if ROI(3)>=im_L
            ROI(3) = im_L;
        end

        ROI(2) = 1 + (j-1)*extractPara.segW + extractPara.thick;
        ROI(4) = ROI(2)+extractPara.segW - 2*extractPara.thick;
        if ROI(4)>=im_W
            ROI(4) = im_W;
        end

        % average neighbor threshold 
        if i<numSegL
            extractPara.trsh = (1/2)*(dynTrsh(i,j) + dynTrsh(i+1,j));
        else
            extractPara.trsh = dynTrsh(i,j);
        end
        % tailor image with ROI
        im_tailor = I(ROI(2):ROI(4),ROI(1):ROI(3));
        rp = extractBubble(im_tailor,ROI,extractPara);
        % collect bubble infomation
        if ~isempty(rp) 
            % area [1], x y [2~3], bounding box [4~7];
            bbinfo_cyc2 = [bbinfo_cyc2;table2array(rp)]; % add to total list
        end
    end
end
% add mask
if ~isempty(bbinfo_cyc2)
    for i=1:length(bbinfo_cyc2(:,1))
        I  = addMask(I,bbinfo_cyc2(i,[2:3,6:7]),extractPara.maskMargin,1); 
    end
end
%% grid detection, thrid cycle, vertical offset
bbinfo_cyc3 = [];    
for i=1:numSegL       
    for j=1:numSegW
        % Segemente the whole picture        
        ROI(1) = 1 + (i-1)*extractPara.segL + extractPara.thick + offsetL;
        ROI(3) = ROI(1) + extractPara.segL - 2*extractPara.thick + offsetL;
        if ROI(3)>=im_L
            ROI(3) = im_L;
        end

        ROI(2) = 1 + (j-1)*extractPara.segW + extractPara.thick + offsetW;
        ROI(4) = ROI(2)+extractPara.segW - 2*extractPara.thick  + offsetW;
        if ROI(4)>=im_W
            ROI(4) = im_W;
        end

        % average neighbor threshold 
        if i<numSegL && j<numSegW
            extractPara.trsh = (1/4)*(dynTrsh(i,j) + dynTrsh(i+1,j)...
                                    + dynTrsh(i,j+1) + dynTrsh(i+1,j+1));
        else
            extractPara.trsh = dynTrsh(i,j);
        end
        % tailor image with ROI
        im_tailor = I(ROI(2):ROI(4),ROI(1):ROI(3));
        rp = extractBubble(im_tailor,ROI,extractPara);
        % collect bubble infomation
        if ~isempty(rp) 
            % area [1], x y [2~3], bounding box [4~7];
            bbinfo_cyc3 = [bbinfo_cyc3;table2array(rp)]; % add to total list
        end
    end
end
if ~isempty(bbinfo_cyc3)
    % add mask
    for i=1:length(bbinfo_cyc3(:,1))
        I  = addMask(I,bbinfo_cyc3(i,[2:3,6:7]),extractPara.maskMargin,1); 
    end
end

% Combine all detected bubbles
bbinfo = [bbinfo;bbinfo_cyc1;bbinfo_cyc2;bbinfo_cyc3];
end