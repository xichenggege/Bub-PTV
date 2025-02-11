function [ind_live] = numLiveBubble(cbbls,frame)
%NUMLIVEBUBBLE returns index of live bubbles 
% That is last frame of bubble is larger than current frame
ind_live = []; i =1;
for kk = 1:numel(cbbls)      % Traverse all recorded bubbles       
    % if bubble is already dead, pass to next bubble
    if (cbbls{kk}.LastFrame >= frame) && (cbbls{kk}.FirstFrame <= frame)
        ind_live(i) = kk;
        i=i+1;
    end  
end
end

