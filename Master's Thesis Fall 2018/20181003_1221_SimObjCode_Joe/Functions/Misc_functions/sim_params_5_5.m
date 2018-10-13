% Waleed Tahir
% 2017-04-04
% Generate parameters to be used in the simulation

% parameters: wavelength and refractive indices
n_medium = 1.33;
n_obj = 1;
lambda=0.6328/n_medium; % (um)

% parameters: da*nzta size
nx = nx;
ny = nx;
nz = 31;
deltaZ = 333; % um (1.54mm)
offsetZ = 26000; % um (17mm)

dpix = 5.5/super_res_factor;  % um
Xlen = nx*dpix; % (110.4 um)
Ylen = nx*dpix; % (110.4 um)


deltaX = Xlen/nx; % um
deltaY = Ylen/ny; % um




