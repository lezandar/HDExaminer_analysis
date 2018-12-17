function color = difference_plot_color_scheme(dat2)
% Specify the ranges of differences in the HX rates and the corresponding
% colors.
% Input: Float or integer. While any real number can be given, only numbers
% between -100 and 100 make biological sense. 
% Output: an array of 3 positive numbers between 0 and 1. Corresponding to
% MATLAB's rgb color set.

if dat2 >=30
     color = [0 0 0.5]; %dark blue
    
elseif dat2 >=15
     color = [0 0.5 1]; % blue
    
elseif dat2 >=5
    color = [0.6 0.93 1];  %light blue
    
elseif dat2 >=-5
    color = [0.6 0.6 0.6]; % gray

elseif (dat2 < -5) && (dat2 > -15)
    color = [1 0.72 0.72];    % pink
    
elseif (dat2 <=-15) && (dat2 >-30)
    color = [1 0 0]; % red
    
elseif (dat2 <=-30)  && (dat2 >-1000)
    color = [0.5 0 0]; % burgundy

else 
    color = [0 1 0]; %green
    
end

