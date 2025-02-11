function [corr] = epipolarCorrelate(bb1,bb2,bb1error,corrPara)
%% EPIPOLARCORRELATE correlate the bubbles of two images
% Start from right image as it has less bubbles

if isempty(bb1) || isempty(bb2)
    corr = [];
else
    for i=1:length(bb2)
        % Inquire points in right image
        epp2 = bb2(i,:);
        % Corresponding epipolar line in left (unregistered image) 
        epiLines = epipolarLine(corrPara.fMatrix',epp2);
        % Search closest point
        D = pdist2(epp2,bb1);
        [dis,index_min] = mink(D,corrPara.trial);
        ind_left = [];

        for j=1:corrPara.trial
            % Inverse transformation from registered to unregistered !!!
            [epp1(1),epp1(2)] = transformPointsInverse(corrPara.tform,...
                bb1(index_min(j),1),bb1(index_min(j),2));

            h = (epp1(1)*epiLines(1)+epiLines(3))/(-epiLines(2)); % height in image 1
            if abs(h-epp1(2))<= bb1error(index_min(j))
                ind_left = index_min(j);
                % delete this point if the distance is small than threshold distance
    %             if dis(j) < corrPara.trshDis
    % %                 bb1(index_min(j),:) = [nan,nan];   
    %             end
                break;
            end
        end
    
        if isempty(ind_left)
            ind_left = nan;
        end   
        corr(i,1) = ind_left;
        corr(i,2) = i;
        corr(i,3) = D(index_min(j));
    end
    % Delete non-matched bubbles
    corr(isnan(corr),:) = [];
end

end




