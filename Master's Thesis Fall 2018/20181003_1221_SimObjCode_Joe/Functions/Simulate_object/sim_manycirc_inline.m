% return object comprising of four slices as mentioned below

function f = sim_manycirc_inline(nx,ny,nz,circ_per_slice)



bndry_offset_px = 20;
circ_size_px = 4;
slice_size_px = nx;

f = zeros(nx,ny,nz);
for i = 1:nz
    f(:,:,1) = circles_slice(circ_per_slice,circ_size_px,slice_size_px,bndry_offset_px);
    f(:,:,3) = f(:,:,1);
    f(:,:,5) = f(:,:,1);
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
