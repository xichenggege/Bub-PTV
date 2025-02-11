function [y_sorted] = sort_data(x,y,x_sort) 

%sort_data sort x with y in the bin of sort_y 
    for i = 1:numel(x_sort)-1
        idx = find(x>x_sort(i) & x<=x_sort(i+1));
        y_sorted{i} = y(idx);
    end
end


