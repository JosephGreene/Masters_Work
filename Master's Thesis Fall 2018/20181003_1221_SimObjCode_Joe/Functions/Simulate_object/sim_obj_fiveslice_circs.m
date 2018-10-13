% return object comprising of four slices as mentioned below

function f = sim_obj_fiveslice_circs(nx,n_obj,n_medium)


bndry_offset_px = 50;
circ_size_px = 12;

slice_size_px = nx;
%% slice1: circles
numcircles = 50;
f(:,:,1) = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px);

%% slice2: circles
numcircles = 45;
f(:,:,2) = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px);

%% slice3: circles
numcircles = 40;
f(:,:,3) = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px);

%% slice4: circles
numcircles = 35;
f(:,:,4) = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px);

%% slice5: circles
f(:,:,5) = zeros(slice_size_px);

%% refractive index value
f = f * (n_obj^2 - n_medium^2);

%% debug
debug = 1;
if(debug)
   figure(1)
   for i = 1:5
      subplot(1,5,i)
      imagesc(f(:,:,i));
      colorbar;
      axis image;
   end
end

end
