% Waleed Tahir
% 2017-04-04
% Generate parameters to be used in the simulation

% keep lateral_n: 256

% parameters: wavelength and refractive indices
n_medium = 1.33;
n_obj = 1;
lambda=0.6328/n_medium; % (um)

% parameters: da*nzta size
nx = lateral_n;
ny = nx;
nz = nz;
deltaZ = deltaZ; % um (1.54mm)
offsetZ = offsetZ; % um (17mm)

dpix = 3.45/8;  % um
Xlen = nx*dpix; % (110.4 um)
Ylen = nx*dpix; % (110.4 um)


deltaX = Xlen/nx; % um
deltaY = Ylen/ny; % um

% this was a test for something, i cant remember what
less_than_pi = (pi*2*deltaX*(Xlen/2))/(lambda*offsetZ);




