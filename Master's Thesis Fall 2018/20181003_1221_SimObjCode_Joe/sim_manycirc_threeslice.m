% return object comprising of four slices as mentioned below

function f = sim_manycirc_threeslice(nx,ny,nz,circ_per_slice,circ_size_px)



bndry_offset_px = 20;
slice_size_px = nx;

f = zeros(nx,ny,nz);
for i = 1:nz
    f(:,:,i) = circles_slice(circ_per_slice,circ_size_px,slice_size_px,bndry_offset_px);
end

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
