function [mask] = Generate_Binary_Axicon(alpha,NAx,NAy,thres)
%{
    alpha adjusts to rate of change of the axicon function
    NAx is the natural aperature grid for the X direction
    NAy is the natural aperature grid for the Y direction
    thres is the binarization thresholding operator
%}

    axicon = exp(-2*pi*1i*alpha*sqrt((NAx.^2+NAy.^2))); %Paper suggests normalizing pupil coordinate
    mask = axicon;
    mask(real(mask) < thres) = -1;
    mask(real(mask) >= thres) = 1;
    mask((NAx.^2+NAy.^2) >= (0.45)^2) = 0;
end