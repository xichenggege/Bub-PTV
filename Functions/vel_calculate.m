function [Traj] = vel_calculate(Traj,step,FPS)
% Velocity calculation with step 
    for i=1:numel(Traj)
        ff = Traj{i}.FirstFrame;
        lf = Traj{i}.LastFrame;
        mm = 1;
        dt = step*1/FPS;  
        
        for kk = (step+ff):1:lf-1
    
           % Note the sequency of position is different from other variables
            ww = kk-ff+1;
    
            Traj{i}.Velocity(mm,:)   = (Traj{i}.pos(ww,:) - Traj{i}.pos(ww-step,:))/dt;     % 3D velocity mm/s
            Traj{i}.Time(mm)         = (1/FPS)*(kk+(1/2)*step);                             % time
            Traj{i}.Cellcenter(mm,:) = (Traj{i}.pos(ww-round(0.5*step),:));  % center of cell (average position)         
            mm = mm + 1;
        end
    end
end
