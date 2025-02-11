function data = stdevFilter(data,stdMutiplier)
%stdevFilter Summary of this function goes here
%##########################################################################
% filter data based on standard deviation
% Delete value beyond ->   mean-n*sigma < data < mean + n*sigma

    mean_val  = mean(data,'omitnan');
    std_val   = std(data,'omitnan');

    lower = mean_val - stdMutiplier*std_val;
    upper = mean_val + stdMutiplier*std_val;

    data(data<lower)=NaN;
    data(data>upper)=NaN;
end

