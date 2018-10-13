% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function f = sim_obj_one(nx,ny,nz,n_obj,n_medium,slicenum)
    f = zeros(nx,ny,nz);
    f(nx/2,ny/2,slicenum) = (n_obj^2 - n_medium^2);
end