function [image] = addMask(image,bbox,mask_margin,maskColorValue)
%ADDMASK Summary of this function goes here
%   Add mask to image 
%   'image' -> image
%   'bbox'  -> bounding box, x,y (center), length and width
% add mask for detected bubbles  
    hL  = bbox(3)/2+mask_margin;
    hW  = bbox(4)/2+mask_margin;
    xc = [bbox(1)-hL, bbox(1)+hL, bbox(1)+hL, bbox(1)-hL];
    yc = [bbox(2)-hW, bbox(2)-hW, bbox(2)+hW, bbox(2)+hW];
    mask  = roipoly(image,xc,yc);
    image(mask) = maskColorValue;  % =1 white in gray image; =255 in normal
end

