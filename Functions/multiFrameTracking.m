function pp = multiFrameTracking(pos,Vlimit,dt)
%multiFrameTracking Summary of this function goes here
%   return predicted point based on 'multi-frames' approach 
%##########################################################################
%--------------------------------------------------------------------------
% 'pos'      -> (x,y) of previous center, n-1, n-2, n-3
% 'Vlimit'   -> limitation of velocity calculation [m/s]
%--------------------------------------------------------------------------
Vlimit = Vlimit*1e3; % m/s -> mm/s

if length(pos(:,1))==1 % if only 1 frame before, use minimum distance
    pp  = pos; % predicted point

elseif length(pos(:,1))==2 % if only 2 frame before, use constant velocity to predict
    % Velocity [m/s]
    vm  = (pos(2,:)-pos(1,:))/dt; 
    if norm(vm)>= Vlimit   % if velocity exceed maximum rising velocity
        vm = vm.*(Vlimit/norm(vm));
    end
    % predicted point
    pp  = pos(2,:) + vm.*dt;

elseif length(pos(:,1))>=3 % 3 points ensures the use of multi-frames prediction
    pos = flipud(pos);   % from start to end -> end to start
    % Velocity [m/s]
    vm1 = (pos(1,:)-pos(2,:))/dt;
    vm2 = (pos(2,:)-pos(3,:))/dt;
    if norm(vm1)>= Vlimit   % if velocity exceed maximum rising velocity
        vm1 = vm1.*(Vlimit/norm(vm1));
    end
    if norm(vm2)>= Vlimit   % if velocity exceed maximum rising velocity
        vm2 = vm2.*(Vlimit/norm(vm2));
    end
    % vm3 = a2*dt + vm2;
    vm0 = 2*vm1-vm2;
    % predicted point
    pp  = pos(1,:) + vm0.*dt;
end
end