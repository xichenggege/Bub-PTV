function [Results] = btv2grid3(Traj,world_time,Lgrid,L,H,D,FPS,fitmethod)
% Distributed bubble informatio into grid
% Each frame has a field
% test

for frame=1:numel(world_time)
    current_time = world_time(frame);     % current moment
    % Create matrix for field data of each frame
    Results{frame}.Uint  = nan*Lgrid;
    Results{frame}.U     = 0*Lgrid;
    Results{frame}.V     = 0*Lgrid;
    Results{frame}.W     = 0*Lgrid;
    Results{frame}.nbb   = 0*Lgrid;       % number of bubbles
    
    %% obtain bubble info in current frame
    xcoord = [];  ycoord = []; zcoord = [];
    u = [];  v = [];  w = [];
    for kk=1:numel(Traj)
        ind_frame = find(abs(Traj{kk}.Time-current_time)<=(1*0.2/FPS)); % introduce a 20% dt as error
        if ~isempty(ind_frame)
            xcoord = [xcoord Traj{kk}.Cellcenter(ind_frame,1)];
            ycoord = [ycoord Traj{kk}.Cellcenter(ind_frame,2)];
            zcoord = [zcoord Traj{kk}.Cellcenter(ind_frame,3)];
            u      = [u     Traj{kk}.Velocity(ind_frame,1)];
            v      = [v     Traj{kk}.Velocity(ind_frame,2)]; 
            w      = [w     Traj{kk}.Velocity(ind_frame,3)]; 
        end
    end
    
    %% Record velocity, both interpolation and non-interpolation
    if numel(xcoord) > 3   % at least 4 points
        % Non-interpolation, record existed bubble only
        for i = 1:numel(xcoord)
            [~,ind_L] = min(abs(L-xcoord(i)));
            [~,ind_H] = min(abs(H-ycoord(i)));
            [~,ind_D] = min(abs(D-zcoord(i)));
            Results{frame}.U(ind_H,ind_L,ind_D)  = Results{frame}.U(ind_H,ind_L,ind_D) + u(i);
            Results{frame}.V(ind_H,ind_L,ind_D)  = Results{frame}.V(ind_H,ind_L,ind_D) + v(i);
            Results{frame}.W(ind_H,ind_L,ind_D)  = Results{frame}.W(ind_H,ind_L,ind_D) + w(i);
            Results{frame}.nbb(ind_H,ind_L,ind_D)= Results{frame}.nbb(ind_H,ind_L,ind_D) + 1;
        end
    end

    if mod(frame,10)==1
    fprintf('======================================= \n');
    fprintf('btv2grid is progressing ... %3.1f%%     \n',frame/numel(world_time)*100);
    fprintf('======================================= \n');
    end 
end

end