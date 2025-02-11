function [cbbls,bblist] = TrackingMatching(I,frame,cbbls,bblist,...
    corrPara,extractPara,trackPara,lastFrame,detector_large)

%==========================================================================
% Tracking and matching
% Tracking previous detected / spatial correlated bubbles with small ROI
%==========================================================================
% Demo version v1.1
% 2022-10-16
% © Xicheng Wang, Yun Feng, Dmitry Grishchenko, Pavel Kudinov
%==========================================================================

%==========================================================================
% 1. Bubble tracking
%==========================================================================
% Demo version v1.0
% 2022-06-17
% © Yun Feng, Xicheng Wang, Dmitry Grishchenko, Pavel Kudinov
%==========================================================================

%% |Step1| -> tracking previous correlated bubbles and remove them from figure
%---------------------------------------------------
% Tracking bubble main progress
%---------------------------------------------------
% find still alive bubbles (lastframe>frame)
ind_live = numLiveBubble(cbbls,frame);
% update bubble information (tracking) after it first frame (when it is observed) 
bbCurrentInfo=[];
% if ~isempty(ind_live)
% Particle tracking with 'Multi frames approach'
for cam=1:2 % left and right camera
    [im_H,im_W] = size(I{cam});
    bbCurrentInfo{cam}=zeros(ind_live(end),7);
    for kk = ind_live      % Traverse all recorded bubbles
        % Extract bubbles from a ROI region
        preframe = length(cbbls{kk}.pos{cam}(:,1));   % frame before current frame
        % Create a 32 x 32 ROI window for local detection
        ROI      = createROIwindow(cbbls{kk}.pos{cam}(preframe,:),16,16,im_H,im_W);

        % tailor image with ROI
        im_tailor = I{cam}(ROI(2):ROI(4),ROI(1):ROI(3));           
        extractPara{cam}.trsh = adapCannyTrsh(im_tailor,0.85,0.95);
        
        % If use RCNN for local Roi detection
        [rp] = localRoiDdetection(im_tailor,ROI,extractPara{cam});  
  
        % kill if not detected
        if isempty(rp)
            cbbls{kk}.LastFrame = frame;
            bbCurrentInfo{cam}(kk,:) = nan;
            continue;
        end 
        
        % Multi frames approach 
        % If multi bubbles in a ROI, how to match them with previous
        % predicted position by last 2 or 3 frames
        pp= multiFrameTracking(cbbls{kk}.pos{cam},trackPara.Vlimit,trackPara.dt);
        % sort nearest point
        [~, index]  = min(pdist2(rp(:,[2:3]),pp));
        bbCurrentInfo{cam}(kk,:) = rp(index,:);      
    end
end

% Update live bubble
ind_live = numLiveBubble(cbbls,frame);

%% |Step2| If bubbles from two images (correlated ROI in previous frame) are stereo correlated 
% Find best possible bubbles in ROI
%==================================================================
% could be further improved
% assume ROI only contains single bubble
%==================================================================
% stereo matching
corr_firstCycle = [];
for kk=ind_live
    epi_error =  (1/2)*norm([bbCurrentInfo{1}(kk,6),bbCurrentInfo{1}(kk,7)]);
    if bubbleSingleCorr(bbCurrentInfo{1}(kk,2:3),bbCurrentInfo{2}(kk,2:3),epi_error,corrPara) % if correlated 
        % save bubble information in current frame
        corr_firstCycle(kk) = 1;
        for cam=1:2
            cbbls{kk}.pos{cam}  = [cbbls{kk}.pos{cam}  ; bbCurrentInfo{cam}(kk,2:3)];
            cbbls{kk}.Area{cam} = [cbbls{kk}.Area{cam} ; bbCurrentInfo{cam}(kk,1)];
            cbbls{kk}.Bbox{cam} = [cbbls{kk}.Bbox{cam} ; bbCurrentInfo{cam}(kk,4:7)];
        end
    else % if not correlated, assume bubble is missing
        cbbls{kk}.LastFrame = frame-1; % not detect in the future 
    end
end
ind_live = numLiveBubble(cbbls,frame);

%% |Step3| -> detection remaining bubbles
%  based on correlated bubbles in current frame
if ~isempty(ind_live)
    for cam = 1:2
        matchedpoint{cam}  = [];       
        for i = 1:numel(ind_live)
            matchedpoint{cam}(i,:) = cbbls{ind_live(i)}.pos{cam}(end,:); 
        end

        bbCurrentInfo{cam}(isnan(bbCurrentInfo{cam}(:,1)),:) = [];
        bbCurrentInfo{cam}(find(bbCurrentInfo{cam}(:,1)==0),:) = [];

        % remove tracked and matched points
        bbNewInfo{cam} = [];
        [~, ia] = setdiff(bbCurrentInfo{cam}(:,[2,3]),matchedpoint{cam},'rows');
        bbNewInfo{cam} = bbCurrentInfo{cam}(ia,:);

        % Add mask
        Imask{cam} = I{cam};
        if ~isempty(bbCurrentInfo{cam})
            for i=1:length(bbCurrentInfo{cam}(:,1))
                Imask{cam} = addMask(Imask{cam},bbCurrentInfo{cam}(i,[2:3,6:7]),extractPara{cam}.maskMargin,1);
            end
        end
        
        % Individual detection using  (grid offset + adaptive threshold)
        % |area| |center X| |center Y| |Bbox: x0 y0 L W|
        [bbNew2,Imask{cam}] =  griddetection(Imask{cam},extractPara{cam}); 

        % For faster RCNN detection (as training input)
        im_detector = uint8(Imask{cam}*255);

        % The remaining part is detected by faster-RCNN
        inputSize = [128 128];
        if ~isempty(detector_large)
            [bbinfo] = FasterRCNNdetection(im_detector,detector_large,inputSize);  
            % Combine grid detection + FasterRCNN detection
            bbNewInfo{cam} = [bbNewInfo{cam};bbNew2;bbinfo]; 
        end   
           
        bbNewInfo{cam} = [bbNewInfo{cam};bbNew2]; 
        bblist{cam}{frame} = [bbCurrentInfo{cam};bbNewInfo{cam}];
    end
end

%% |Step4| -> Stereo correlated remaining bubbles and add as new
% Local based error for epipolar constraint
epi_error = (1/2)*sqrt(bbNewInfo{1}(:,6).^2 + bbNewInfo{1}(:,7).^2);

bb1_Reg = [];
bb2_Reg = [];
% Transform bubble to registered coordinate 
[bb1_Reg(:,1),bb1_Reg(:,2)] = transformPointsForward(corrPara.tform,bbNewInfo{1}(:,2),bbNewInfo{1}(:,3));
bb2_Reg = bbNewInfo{2}(:,2:3); % Only left image is registered 

% Epipolar constraint
corrEpipolar  = epipolarCorrelate(bb1_Reg,bb2_Reg,epi_error,corrPara);
if isempty(corrEpipolar)
else
    % Add new bubbles to 'cbbls'
    listend = numel(cbbls);
    for kk=1:numel(corrEpipolar(:,1)) 
        for cam=1:2 % left and right camera
            cbbls{kk+listend}.pos{cam}  = bbNewInfo{cam}(corrEpipolar(kk,cam),2:3);
            cbbls{kk+listend}.Area{cam} = bbNewInfo{cam}(corrEpipolar(kk,cam),1);
            cbbls{kk+listend}.Bbox{cam} = bbNewInfo{cam}(corrEpipolar(kk,cam),4:7);
        end
        cbbls{kk+listend}.listdata     = kk+listend;    % n th bubble
        cbbls{kk+listend}.FirstFrame   = frame;
        cbbls{kk+listend}.LastFrame    = lastFrame;
    end
end

%% Delete bubbles with too short life
frameLength = [];
for kk=1:numel(cbbls)
    frameLength(kk) = cbbls{kk}.LastFrame-cbbls{kk}.FirstFrame;
end

ind_shortLife = find(frameLength<=trackPara.limitFrameL);
cbbls(ind_shortLife) = [];                 % avoid out of memory, delete too short bubbles

end
    


