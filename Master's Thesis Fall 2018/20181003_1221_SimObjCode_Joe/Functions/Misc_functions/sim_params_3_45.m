% Waleed Tahir
% 2017-04-04
% Generate parameters to be used in the simulation

% parameters: wavelength and refractive indices
n_medium = 1.33;
n_obj = 1;
lambda=0.6328/n_medium; % (um)

% parameters: da*nzta size
nx = nx*super_res_factor;
ny = nx;
nz = 5;
deltaZ = 6000; % um (1.54mm)
offsetZ = 14000; % um (17mm)

dpix = 3.45/super_res_factor;  % um
Xlen = nx*dpix; % (110.4 um)
Ylen = nx*dpix; % (110.4 um)


deltaX = Xlen/nx; % um
deltaY = Ylen/ny; % um

% this was a test for something, i cant remember what
less_than_pi = (pi*2*deltaX*(Xlen/2))/(lambda*offsetZ);




