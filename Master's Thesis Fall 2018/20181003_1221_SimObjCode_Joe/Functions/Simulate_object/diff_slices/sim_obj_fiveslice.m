% return object comprising of four slices as mentioned below

function f = sim_obj_fiveslice(nx,n_obj,n_medium)

slice_size_px = nx;

%% slice1: circles
numcircles = 30;
circ_size_px = 1;
bndry_offset_px = 5;
f(:,:,1) = circles_slice(numcircles,circ_size_px,slice_size_px,bndry_offset_px);

%% slice2: text
xpos = 70;
ypos = 150;
f(:,:,2) = text_slice('PI',slice_size_px,xpos,ypos);

%% slice 3: lines
hor_xpos = [50 100];
hor_ypos = [50 200];
hor_len = [50 75];
ver_xpos = [200];
ver_ypos = [50];
ver_len = [50];
linwidth_px = 1;

f(:,:,3) = lines_slice(hor_xpos,hor_ypos,hor_len,ver_xpos,ver_ypos,...
    ver_len,slice_size_px,linwidth_px);

%% slice4: Xs
xpos = [80 176 176];
ypos = [80 80 176];
f(:,:,4) = text_slice('X',slice_size_px,xpos,ypos);

%% slice 5: empty
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
