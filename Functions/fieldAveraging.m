function [meddata avgdata,stddata,u_avgdata] = fieldAveraging(field,filterMethod,stdMultiplier)
%#############################################################################
% Post processing of field data
% performe filtering 'stdFilter' to remove abnormal data
% Get median value of each grid point along all the frame
% Get averaged value of each grid point along all the frame

for f=1:numel(field) 
    data(:,f) = field{f}(:);       % converte to 1D vector
end

% use filter to remove abonrmal data
if strcmp(filterMethod,'stdFilter')  %% more fiter could be introduced
    for i = 1:length(data(:,1))
        data(i,:)  = stdevFilter(data(i,:),stdMultiplier);
        Neff(i,1)  = numel(find(~isnan(data(i,:))));
    end
else
    fprintf('warnning: no filtering, please check');
end

% calculate median and average data
meddata  = median(data,2,'omitnan'); % '2' -> column wise vector median
avgdata  = mean(data,2,'omitnan');   % '2' -> column wise vector mean
stddata  = std(data,[],2,'omitnan');    % '2' -> column wise standard deviation 

u_avgdata = stddata./sqrt(Neff);


% Reshape the matrix
meddata   =  reshape(meddata,size(field{1}));
avgdata   =  reshape(avgdata,size(field{1}));
stddata   =  reshape(stddata,size(field{1}));
u_avgdata =  reshape(u_avgdata,size(field{1}));
end