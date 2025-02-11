function [meanimg] = preProcessing(raw_images,NumFrames,StereoCameraParam,submean)
%% Subtract mean of all images to remove the background of each image
% It improves results of PTV

meanimg{1}=zeros(size(raw_images{1}{1},1),size(raw_images{1}{1},2)); % initialization
meanimg{2}=zeros(size(raw_images{2}{1},1),size(raw_images{2}{1},2)); % initialization

% Subtract mean of all images to remove the background of each image
if submean
    for cam=1:2
        for frame=1:NumFrames 
            % read and undistort image
            if cam==1
                image = undistortImage(raw_images{1}{frame},StereoCameraParam.CameraParameters1);   % un-distorted
            else
                image = undistortImage(raw_images{2}{frame},StereoCameraParam.CameraParameters2);   % un-distorted
            end
    
            meanimg{cam}=(((frame-1).*meanimg{cam})+double(image))./frame;
    
            if mod(frame,10)==1
            fprintf('Camera %1.0f is progressing ... %3.1f%% \n',cam,(frame-1)/NumFrames*100);
            end
        end     
    end
end

end

