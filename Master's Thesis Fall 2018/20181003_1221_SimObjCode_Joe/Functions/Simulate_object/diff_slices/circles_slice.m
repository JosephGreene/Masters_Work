% this function puts numcircles no. of circles in random positions in
% a 2D slice and returns the slice
function slice = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px)

% generate random indexes for circles
idx_range = [bndry_offset_px,slice_size_px-bndry_offset_px];
idx_x = randi(idx_range,[1 numcircles]);
idx_y = randi(idx_range,[1 numcircles]);
% remove indexes near boundary
idx_x(idx_x<2*circ_size_px) = 0;
idx_y(idx_y<2*circ_size_px) = 0;
% make indexes negative to positive instead of one to positive
idx_x = idx_x - slice_size_px/2;
idx_y = idx_y - slice_size_px/2;

size = slice_size_px;
diameter = circ_size_px/size;

value = ones(numcircles,1);
vertlen = diameter * ones(numcircles,1);
horzlen = diameter * ones(numcircles,1);
xcor = idx_x';
ycor = idx_y';
phi = zeros(numcircles,1);

E = [value vertlen horzlen (xcor/size)*2 (ycor/size)*2 phi];

slice = phantom(E, size);
slice(slice>1) = 1;

end