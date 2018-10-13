% return object comprising of four slices as mentioned below

function f = sim_circ_two(nx,ny,nz,n_obj,n_medium,~)


bndry_offset_px = 50;
circ_size_px = 12;
slice_size_px = nx;
numcircles = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slice 1
% generate random indexes for circles
idx_range = [bndry_offset_px,slice_size_px-bndry_offset_px];
idx_x = randi(idx_range,[1 numcircles]);
idx_y = randi(idx_range,[1 numcircles]);
idx_x = nx/2; %center
idx_y = ny/2; %center
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

f = zeros(nx,ny,nz);
f(:,:,1) = phantom(E, size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slice 2
% generate random indexes for circles
idx_range = [bndry_offset_px,slice_size_px-bndry_offset_px];
idx_x = randi(idx_range,[1 numcircles]);
idx_y = randi(idx_range,[1 numcircles]);
idx_x = 3*nx/4; %center
idx_y = 3*ny/4; %center
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

f(:,:,2) = phantom(E, size);



%% refractive index value
f = f * (n_obj^2 - n_medium^2);

%% debug
debug = 0;
if(debug)
   figure(1)
   for i = 1:nz
      subplot(1,5,i)
      imagesc(f(:,:,i));
      colorbar;
      axis image;
   end
end

end
