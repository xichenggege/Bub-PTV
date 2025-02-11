function [trsh] = adapCannyTrsh(image,multiplier,upperlimit)
%adapCannyTrsh Summary of this function goes here
%   Use canny edge detector to determine the adapative threshold value
%##########################################################################
% 'image'      -> image
% 'multiplier' -> control the mean value
% 'upperlimit' -> upper limit of threshold value
% local adapative thershold value based on 'cannoy' edge detector

trsh = image(edge(image,'canny'));
trsh(trsh==1) = [];   % exclude 'white mask'
trsh = multiplier*mean(trsh);     % threshold multiplier 

% aovid too high threshold value to detect non-exsit bubble (white)
if trsh >=upperlimit
        trsh = upperlimit;
end

end