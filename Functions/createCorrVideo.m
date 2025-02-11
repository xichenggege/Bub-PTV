function [] = createCorrVideo(raw_images,cbbls,videoSetup,StereoCameraParam)
%createCorrVideo creates video to display the correlation of bubbles
%   Detailed explanation goes here
%   Display correlated reuslts in all frame
% Create video
vdtemp = VideoWriter('Videos\tracking_matching');
vdtemp.FrameRate = videoSetup.frameRate;
open(vdtemp);
set(0,'DefaultFigureVisible', 'off');

for frame = videoSetup.start:videoSetup.frameSkip:videoSetup.end
    % Read figures
    im1_distort = raw_images{1}{frame}; 
    im2_distort = raw_images{2}{frame};    
    image{1} = undistortImage(im1_distort,StereoCameraParam.CameraParameters1);   % un-distorted
    image{2} = undistortImage(im2_distort,StereoCameraParam.CameraParameters2);
    % Plot start

    FigHandle = figure; 
    set(FigHandle, 'Position', [50, 50, 2400, 700]);  % keep the same frame size
    for cam=1:2
        % Show left figure
        pos = [0.05+0.5*(cam-1) 0.05 0.45 0.9];
        subplot('Position',pos); hold on
        imshow(image{cam},'Border','tight');   
        % Bubble centers
        ind_live = numLiveBubble(cbbls,frame);
        matchedpoint{cam} = [];
        for i = 1:numel(ind_live)
            currentFrame = frame - cbbls{ind_live(i)}.FirstFrame + 1;
            matchedpoint{cam}(i,:) = cbbls{ind_live(i)}.pos{cam}(currentFrame,:); 
        end
        plot(matchedpoint{cam}(:,1),matchedpoint{cam}(:,2),'+g'); 
    end
    % title
    sgtitle(strcat('bubble detection in frame=',num2str(frame)));

    % Record figure
    F= getframe(gcf);
    writeVideo(vdtemp,F);
    close(FigHandle);

    % Processing bar
    % print output information
    if mod(frame,10)==1
    fprintf('======================================= \n');
    fprintf('Creating video is progressing ... %3.1f%% \n',(frame)/videoSetup.end*100);
    fprintf('======================================= \n');
    end
end
    close(vdtemp);
end