%##########################################################################
% Input File
% Prepare all required parameters in this file
%##########################################################################
% Demo version v1.0
% 2022-06-14
% © Yun Feng & Xicheng Wang 
%==========================================================================
%% If data is available
corrPara_available  = true;
bblist_available    = false;
cbbls_available     = false;
Traj3D_available    = false;
griddata_available  = false;
creatVideo          = true;

%% If subtract mean image
submean = true;


%% Calibration parameters based on 'Camera_Calibration.m'
dataLoc=strcat(num2str(pwd));
% 'StereoCameraParam.mat' from stereo camera calibration
% Refer to 'https://se.mathworks.com/help/releases/R2024b/vision/camera-calibration.html?searchHighlight=stereo+camera+calibration&s_tid=doc_srchtitle'
file = load(strcat(dataLoc(1:end),'\Data\pre-processing\StereoCameraParam.mat'));   
StereoCameraParam = file.camParam;
RefPlane          = file.RefPlane; 

%% Video information
% This is a demo code
% 200 frames are used which were extracted from original videos
dataLoc = strcat(num2str(pwd));
% Read all images, overall 200 frames

file          = load(strcat(dataLoc,'\Data\pre-processing\videoImages_1.mat'));
imageFile{1} = file.im';
file         = load(strcat(dataLoc,'\Data\pre-processing\videoImages_2.mat'));
imageFile{2} = file.im';
file         = load(strcat(dataLoc,'\Data\pre-processing\videoImages_3.mat'));
imageFile{3} = file.im';
file         = load(strcat(dataLoc,'\Data\pre-processing\videoImages_4.mat'));
imageFile{4} = file.im';
raw_images   = [];
% Combine them together
k = 1;
for i = 1:4
    file  = load(sprintf('%s\\Data\\pre-processing\\videoImages_%d.mat',dataLoc,i));
    for j = 1:numel(file.im{1})
        raw_images{1}{k} = file.im{1}{j};
        raw_images{2}{k} = file.im{2}{j};
        k = k+1;
    end
end

% frame range
NumFrames    = k-1;   % number of total frames
lastFrame    = NumFrames;          
startFrame   = 1;

%% Bubble Tracking

% Region of Interest (ROI),
ROIHalfW   = 8;   % half width of ROI window
ROIHalfH   = 8;   % half length of ROI window 

% frame rate
FPS = 6300;

% for stagation test, this is maximum possible velocity [m/s]
trackPara.Vlimit       = 5;  % use effective velocity to limit multi-frame tracking prediction
trackPara.dt           = 1/FPS; 
trackPara.limitFrameL  = 10; % Ignore short frame length bubble
trackPara.searchRadius = 2; % nearest searching radius [pixel]
%% Stereo matching
% Bubble correlation or matching
corrPara.trial   = 2;    % maximum number of trial of searching

%% Bubble detecting/extracting
% Subtract mean of all images to remove the background of each image
submean = true;

% Left camera
extractPara{1}.trsh=0.2;           % initial thershhold value
extractPara{1}.trsh_max = 0.95;    % maximum threshold value
extractPara{1}.areafilt= [6 65];  % area filter
extractPara{1}.segL    = 6*ROIHalfH;   % segement length 
extractPara{1}.segW    = 6*ROIHalfW;   % segement Witdth 
extractPara{1}.thick   = 3;     % thickness of grid mask
extractPara{1}.maskMargin = 1;  % mask margin
% Right camera
extractPara{2}.trsh=0.2;     
extractPara{2}.trsh_max = 0.95;    % maximum threshold value
extractPara{2}.areafilt= [6 65];
extractPara{2}.segL    = 6*ROIHalfH;   % segement length 
extractPara{2}.segW    = 6*ROIHalfW;   % segement Witdth 
extractPara{2}.thick   = 3;        % thickness of grid mask
extractPara{2}.maskMargin = 1;  % mask margin

% Load fasterRCNN detector
% No faster RCNN for water test
% detector_large = importdata('Data\faster_rcnn_checkpoint__7276__2022_12_09__20_18_57.mat');
detector_large = [];

%% PostProcessing
% Velocity calculation
tlag         = 10;      % velocity calculation per 'timelag' intervals (better to be odd number)
world_time   = [1:1:lastFrame]*1/FPS;

