%##########################################################################
% Bubble based Particle Tracking Velocimetry (Bub-PTV) 
%##########################################################################
%==========================================================================
% To perform analysis of stereo videos as:
% 0. Read input file
% 1. Bubble detection
% 2. Tracking and stereo matching
% 3. Post-processing (3D trajectory, velocity and btv2grid etc)
% 4. Further post-processing by user defined
%==========================================================================
% Demo version v1.3
% 2022-10-16
% © Xicheng Wang, Yun Feng, Dmitry Grishchenko, Pavel Kudinov
%==========================================================================

clc
clear all

%% Load input file
% Load video & calibration matrix of camera
InputFileName = char('InitialSetup_SEFW7_V1');
try
    run(InputFileName);
catch
    disp(strcat('Error:',InputFileName,...
        'is missing, please put it under the same path'));
end

% add main code path to include all functions
addpath('.\Functions');

figurePre_path = strcat(pwd,'\Figures\pre-processing\');
figurePro_path = strcat(pwd,'\Figures\processing\');
mkdir(figurePre_path);
mkdir(figurePro_path);

dataFolder = strcat('Data\post-processing\');
mkdir(dataFolder);

%----------------------------------------------------------------------
%% Tracking main code
%----------------------------------------------------------------------

if ~cbbls_available
    %----------------------------------------------------------------------
    %% Initialization
    %----------------------------------------------------------------------
    [meanimg] = preProcessing(raw_images,NumFrames,StereoCameraParam,submean);
    if submean
        imwrite(mat2gray(meanimg{1}),[figurePre_path '\Mean-cam1.png'],'png');
        imwrite(mat2gray(meanimg{2}),[figurePre_path '\Mean-cam2.png'],'png');
    end

    %----------------------------------------------------------------------
    %% Initial matching
    %----------------------------------------------------------------------
    frame = startFrame;
    for cam = 1:2   
        % load image
        if cam == 1
            im_distort = raw_images{1}{frame};     % load  
            I{cam}       = undistortImage(im_distort,StereoCameraParam.CameraParameters1);   % un-distorted                         
        elseif cam == 2
            im_distort = raw_images{2}{frame};     % load  
            I{cam}       = undistortImage(im_distort,StereoCameraParam.CameraParameters2);   % un-distorted
        end
        % global detection 
        [bblist{cam}{frame}] = global_detection(I{cam},submean,...
                meanimg{cam},extractPara{cam},detector_large); 
    end
    
    [cbbls,corrPara] = InitialMatching(I,bblist,StereoCameraParam.FundamentalMatrix,...
        startFrame,lastFrame,corrPara_available);
    
    %----------------------------------------------------------------------
    %% Tracking and matching
    %----------------------------------------------------------------------
    for frame = startFrame+1:1:lastFrame
        for cam = 1:2   
            % load image
            if cam == 1
                im_distort = raw_images{1}{frame};     % load  
                I{cam}       = undistortImage(im_distort,StereoCameraParam.CameraParameters1);   % un-distorted                         
            elseif cam == 2
                im_distort = raw_images{2}{frame};     % load  
                I{cam}       = undistortImage(im_distort,StereoCameraParam.CameraParameters2);   % un-distorted
            end

            % Subtract mean of all images to remove the background of each image
            %It improves results of SBTV
            if submean == 1 
                I{cam} = mat2gray(double(I{cam})-meanimg{cam});  % subtract mean value to remove noisy background
            else
                I{cam} = mat2gray(I{cam});  % subtract mean value to remove noisy background
            end
        end
      
        [cbbls,bblist] = TrackingMatching(I,frame,cbbls,bblist,...
        corrPara,extractPara,trackPara,lastFrame,detector_large);
    
        % Processing bar
        % print output information
        if mod(frame,10)==1
            fprintf('======================================= \n');
            fprintf('Tracking is progressing ... %3.1f%% \n',(frame)/lastFrame*100);
            fprintf('======================================= \n');
            
            for cam = 1:2
                bBoxPlot = insertShape((I{cam}),'rectangle',bblist{cam}{frame}(:,4:7), 'LineWidth', 1, 'Color','blue');
                name=sprintf('%scam%d_BubbleDetection_frame%d.png',figurePro_path,cam,frame);
                imwrite(bBoxPlot,name);
            end

        end
    end
    % Save data
    save(strcat(dataFolder,'bblist'),'bblist'); % Bubbles separately detected from two cameras
    save(strcat(dataFolder,'cbbls'),'cbbls');   % Correlated & tracked bubbles
   
else
    % Load directly
    file    = load(strcat(dataFolder,'cbbls')); % load it
    cbbls   = file.cbbls;
end

%----------------------------------------------------------------------
%% create video for preview your tracking & matching 
%----------------------------------------------------------------------
if creatVideo
    videoSetup.frameRate = 5;
    videoSetup.frameSkip = 1;
    videoSetup.start     = startFrame;
    videoSetup.end       = 200;
    createCorrVideo(raw_images,cbbls,videoSetup,StereoCameraParam);
end

%----------------------------------------------------------------------
%%  Calculate 3D trajectory of all matched bubbles
%----------------------------------------------------------------------
smooth_window = 10;
if ~Traj3D_available
    k=1;  
    for i=1:1:numel(cbbls)   
        % create 3D trajectory 
        [Traj{i}.pos Traj{i}.FirstFrame Traj{i}.LastFrame] = ...
            compute3DTrajetory(cbbls{i},StereoCameraParam,RefPlane,'rotate',smooth_window);
        % Processing bar
        % print output information
        if mod(i,10)==1
        fprintf('======================================= \n');
        fprintf('Trajetories is progressing ... %3.1f%%  \n',(i)/numel(cbbls)*100);
        fprintf('======================================= \n');
        end 
    end  

    Hoff= 0; % injection hole vertical offset, estimated by eyes
    Loff= 0;
    Doff= 0;
    for kk=1:numel(Traj)
        Traj{kk}.pos(:,2) = -Traj{kk}.pos(:,2);  % Flip vertical coordinate !
    end   
    save(strcat(dataFolder,'Traj'),'Traj');   % Correlated & tracked bubbles
else
    file    = load('Data\post-processing\Traj3'); % load it
    Traj    = file.Traj;
end

