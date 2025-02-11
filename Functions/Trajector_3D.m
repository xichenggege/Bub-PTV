%##########################################################################
% Calculation of world trajetory of two bubbles
% 1. 'Rotate' means rotate the trajeotry based on reference plane
%##########################################################################
%==========================================================================
% 2022-02-02
% © Yun Feng & Xicheng Wang 
%==========================================================================
function [traj,ff,lf] = Trajector_3D(bb1,bb2,StereoCameraParam,text1)

    smooth_traj=1;         % if smooth the trajetory
    smooth_window = 100;   % span of smoothing window
    % Load stereoCamera parameters
    RefPlane = StereoCameraParam.RefPlane;
    camParam = StereoCameraParam.camParam;

    % Find real last frame
    ff1 = bb1.FirstFrame;
    lf1 = bb1.LastFrame;

    ff2 = bb2.FirstFrame;
    lf2 = bb2.LastFrame;

    ff = max(ff1,ff2);
    lf = min(lf1,lf2);

    X1  = bb1.pos(ff:lf,1);
    Y1  = bb1.pos(ff:lf,2);
    X2  = bb2.pos(ff:lf,1);
    Y2  = bb2.pos(ff:lf,2);

    % Smooth the trajectories
    if smooth_traj % if smooth

        X1 = smooth(X1, smooth_window,'rloess');
        Y1 = smooth(Y1, smooth_window,'rloess');
        
        X2 = smooth(X2, smooth_window,'rloess');
        Y2 = smooth(Y2, smooth_window,'rloess');
    end
   
    % Undistort the points according to information from calibration
    XY1  = undistortPoints([X1 Y1], camParam.CameraParameters1);
    XY2  = undistortPoints([X2 Y2], camParam.CameraParameters2); 
    traj_raw = triangulate(XY1, XY2, camParam);    

%     plot3(traj_raw(:,1),traj_raw(:,2),traj_raw(:,3)); hold on
%     plot3(traj(:,1),traj(:,2),traj(:,3)); hold on

    if isequal(text1,'rotate')  % isequal(a,b) return 1 if a=b
        x =   RefPlane(:,1);
        y =   RefPlane(:,2);
        z =   RefPlane(:,3);
        
        X = [ones(numel(x),1), x];
        by = (X'*X)^-1*X'*z;
        Results.alphaY = atan(by(2)); % rad         
        
        Y = [ones(numel(y),1), y];
        bx = (Y'*Y)^-1*Y'*z;                    
        Results.alphaX = atan(bx(2)); % rad 
        
        bz = (X'*X)^-1*X'*y;                    
        Results.alphaZ = atan(bz(2));   

        alphaY =  Results.alphaY*(abs( Results.alphaY)>2*pi/180);
        alphaX =  Results.alphaX*(abs( Results.alphaX)>2*pi/180);
        alphaZ =  Results.alphaZ*(abs( Results.alphaZ)>2*pi/180);  
        % rotation around Z axis
        AZ = [cos(alphaZ) -sin(alphaZ) 0;...
             sin(alphaZ)  cos(alphaZ) 0;...
            0               0             1];
        
        % rotation around Y axis
        AX = [cos(alphaY) 0 sin(alphaY);...
            0              1 0          
            -sin(alphaY) 0 cos(alphaY)];        
        
        % rotation around X axis
        AY = [1              0            0
             0              cos(alphaX) -sin(alphaX);...
             0              sin(alphaX)  cos(alphaX)];         
         
        % rotation matrix A 
        A = AX*AY*AZ;    
        traj = (A*traj_raw')';
    end
end