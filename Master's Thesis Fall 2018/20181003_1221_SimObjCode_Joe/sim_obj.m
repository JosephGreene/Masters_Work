%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate datacube for holography
% Waleed Tahir
% Date: 3 Jan 2018
function f = sim_obj(nx,ny,nz,par_radius_px,circ_per_slice_vec)
    %add_paths;

%     nx = 401;
%     ny = 401;
%     nz = 10;
%     par_radius_px = 2;
%     circ_per_slice_vec = [10 50];

%     f = zeros(nx,ny,nz);

    for circ_per_slice = circ_per_slice_vec

        f = sim_manycirc_threeslice(nx,ny,nz,circ_per_slice,par_radius_px);

%         fn = sprintf('gtruth/numperslice_%d/slices',circ_per_slice);
%         mkdir(fn);
%         for i = 1:nz
%             h = imagesc(f(:,:,i));
%             axis image;
%             filename = sprintf('%s/slice%d.png',fn,i);
%             saveas(h,filename,'png');
%         end 

        f_obj = f;
        %save(sprintf('gtruth/numperslice_%d/f_obj_%d',circ_per_slice,...
            %circ_per_slice),'f_obj');
    end
end 