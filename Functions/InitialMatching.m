function [cbbls,corrPara] = InitialMatching(I,bblist,fMatrix,startFrame,lastFrame,corrPara_available)
%% Image registration
unregistered = I{1};    % register left image
ortho        = I{2};

%% Find at least 8 matched points
if corrPara_available
    try
        % Parameters for correlation, obtained by manually selected at least 8 matched points
        % Details refer to "https://se.mathworks.com/help/releases/R2024b/images/ref/cpselect.html?searchHighlight=cpselect&s_tid=doc_srchtitle"
        file = load('Data\pre-processing\corrPara');
    catch
        disp(strcat('Error: "corrPara" is missing, please put it under \Data\post-processing'));
    end
    corrPara = file.corrPara;
    mp1 = corrPara.mp1;
    mp2 = corrPara.mp2;
else % Can switch 'corrPara_available=false' in 'InitialSetup_SEFW7_V1.m'
    mp1 = zeros(1,2);
    while numel(mp1(:,1))<8
        [mp1,mp2] = cpselect(unregistered,ortho,'Wait',true);
        if numel(mp1(:,1))<=8
            display('Please select at least 8 matched points');
        end
    end
end

%% 'projective' transformation type    
tform = fitgeotrans(mp1,mp2,'projective'); 
Rfixed = imref2d(size(ortho));
registered = imwarp(unregistered,tform,'OutputView',Rfixed);

% Unregistered image
figure
set(gcf, 'Position', [50, 300, 600, 510]);
imshowpair(unregistered,ortho,'falsecolor'); hold on
title('Image unregistration');
% Image show for test
figure
set(gcf, 'Position', [50, 300, 600, 510]);
imshowpair(registered,ortho,'falsecolor'); hold on
title('Image registration');

%% Epipolar line on registered image
% Directly from Stereo Calibration matrix (not used in this test)
% Epipolar line for display

% Compute fundemental matrix !!!
% inliers returns index of points used for calculation (valid) 
% at least 8 valid points are required

% % Transform bubble centers from raw images
% [mp1_Reg(:,1),mp1_Reg(:,2)] = transformPointsForward(tform,mp1(:,1),mp1(:,2));
% mp2_Reg = mp2; % Only image 1 is registered 
% [fMatrix, inliers] = estimateFundamentalMatrix(mp1_Reg,...
%     mp2_Reg,'Method','Norm8Point'); 

epiLines = epipolarLine(fMatrix',mp2);
points   = lineToBorderPoints(epiLines,size(unregistered));

%% Correlation according to epipolar constraint
% Add bubbles from unregistered images (first frame)
bbinfo{1} = bblist{1}{startFrame};
bbinfo{2} = bblist{2}{startFrame};
% Transform bubble centers from raw images
[bb1_Reg(:,1),bb1_Reg(:,2)] = transformPointsForward(tform,bbinfo{1}(:,2),bbinfo{1}(:,3));
bb2_Reg = bbinfo{2}(:,2:3); % Only left image is registered 
% epipolar constraint, 0.5*diagonal length
ep_error = 0.5*sqrt(bbinfo{1}(:,6).^2+bbinfo{1}(:,7).^2);

% Correlation parameters
corrPara.fMatrix    = fMatrix; % fundamental matrix
corrPara.tform      = tform;   % transformation matrix of image registration
corrPara.mp1        = mp1;
corrPara.mp2        = mp2;
corrPara.trial      = 2;    % maximum number of trial of searching
% save('Data\corrPara','corrPara'); % Save 'Correlation parameters', next time won't do this again

% Epipolar constraint
corrEpipolar  = epipolarCorrelate(bb1_Reg,bb2_Reg,ep_error,corrPara);
% % distance for absolute correct 
% corrPara.trshDis  = mean(mink(corrEpipolar(:,3),round(numel(corrEpipolar(:,3))*0.5)));

% Saved correlated bubbles
for kk=1:numel(corrEpipolar(:,1)) 
    for cam=1:2 % left and right camera
        cbbls{kk}.pos{cam}  = bbinfo{cam}(corrEpipolar(kk,cam),2:3);
        cbbls{kk}.Area{cam} = bbinfo{cam}(corrEpipolar(kk,cam),1);
        cbbls{kk}.Bbox{cam} = bbinfo{cam}(corrEpipolar(kk,cam),4:7);
    end
    cbbls{kk}.FirstFrame  = startFrame;
    cbbls{kk}.LastFrame   = lastFrame;
    cbbls{kk}.listdata    = kk;
end

% For display
figure
showMatchedFeatures(registered,ortho,...
    bb1_Reg(corrEpipolar(:,1),:),bb2_Reg(corrEpipolar(:,2),:),'falsecolor');
title('Original images and matching feature points by registration and epipoline');

% print output information
fprintf('%3.0f pairs are correlated with a total number of %3.0f \n',numel(corrEpipolar(:,1)),length(bbinfo{2}(:,1)));
end

