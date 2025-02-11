function [Results] = btv2grid2(Traj,world_time,lgrid,L,R,FPS)
% Distributed bubble informatio into grid
% Each frame has a field
% test

for frame=1:numel(world_time)
    current_time = world_time(frame);     % current moment 
    Results{frame}.U_avg2= 0*lgrid;
    Results{frame}.V_avg2= 0*lgrid;
    Results{frame}.nbb2  = 0*lgrid;

    %% obtain bubble info in current frame
    xcoord = [];  ycoord = []; zcoord = [];
    u = [];  v = [];  w = [];
    for kk=1:numel(Traj)
        ind_frame = find(abs(Traj{kk}.Time-current_time)<=(1*0.2/FPS)); % introduce a 20% dt as error
        if ~isempty(ind_frame)
            xcoord = [xcoord Traj{kk}.Cellcenter(ind_frame,1)];
            ycoord = [ycoord Traj{kk}.Cellcenter(ind_frame,2)];
            zcoord = [zcoord Traj{kk}.Cellcenter(ind_frame,3)];
            u      = [u     Traj{kk}.Velocity(ind_frame,1)];   % x velocity
            v      = [v     Traj{kk}.Velocity(ind_frame,2)];   % y velocity
        end
    end
    

    % Record velocity, both interpolation and non-interpolation
    if numel(xcoord) > 3   % at least 4 points
        % Non-interpolation, record existed bubble only
        % radial distance
        rcoord = (ycoord./abs(ycoord)).*sqrt(zcoord.^2+ycoord.^2);
        for i = 1:numel(xcoord)
            [~,ind_L] = min(abs(L-xcoord(i)));
            [~,ind_R] = min(abs(R-rcoord(i)));
    
            Results{frame}.U_avg2(ind_R,ind_L)  = Results{frame}.U_avg2(ind_R,ind_L) + u(i);
            Results{frame}.V_avg2(ind_R,ind_L)  = Results{frame}.V_avg2(ind_R,ind_L) + v(i);
            Results{frame}.nbb2(ind_R,ind_L)    = Results{frame}.nbb2(ind_R,ind_L) + 1;
        end
    end


    if mod(frame,10)==1
    fprintf('======================================= \n');
    fprintf('btv2grid2 is progressing ... %3.1f%%     \n',frame/numel(world_time)*100);
    fprintf('======================================= \n');
    end 
end

end