% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function f = sim_obj_center(nx,ny,nz,n_obj,n_medium,num_particles)


    %% create num_particles # of bubbles in a 256x256x20 cube
    bndry = 25;
    f_cent = zeros(nx-(2*bndry),ny-(2*bndry),nz);
    [m,n,~] = size(f_cent);
    bubble_idx = randi(m*n*nz,[1 num_particles]);
    f_cent(bubble_idx) = n_obj^2 - n_medium^2;
    
    % put the 256x256x20 cube above in the middle of the 1024x1024x20 cube
    f = zeros(nx,ny,nz);
    f(bndry:bndry+(m-1),bndry:bndry+(n-1),:) = f_cent;
    
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