%##########################################################################
% Calculation of world trajetory of two bubbles
% Applied to 'SEF_ImPro_v8'
% 1. 'Rotate' means rotate the trajeotry based on reference plane
%##########################################################################
%==========================================================================
% 2022-05-16
% © Yun Feng & Xicheng Wang 
%==========================================================================
function [traj,ff,lf] = compute3DTrajetory(bbl,StereoCameraParam,RefPlane,text1,smooth_window)
   
    % find avaliable data
    ff = bbl.FirstFrame;
    lf = bbl.LastFrame;
    for cam=1:2
        X{cam} = smooth(bbl.pos{cam}(:,1), smooth_window,'rloess');
        Y{cam} = smooth(bbl.pos{cam}(:,2), smooth_window,'rloess');
        XY{cam} = [X{cam},Y{cam}];
%         % not necessary as the image has been undistorted 
%         XY_undistort{cam}  = undistortPoints([X{cam},Y{cam}], StereoCameraParam.CameraParameters1);
    end

    traj = triangulate(XY{1}, XY{2},StereoCameraParam);  

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
        traj = (A*traj')';
    end
end