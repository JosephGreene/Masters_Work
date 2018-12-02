F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)
pre = @(x) (x)/norm(x); %Preprocessing step to normalize input data for autocorrelation


%Physical parameters of miniscope (unchagable) in real space
NA0 = 0.45; %NA of the system -> range of ATF (2*NA is range of OTF)
f = 3.3; %mm
lambda = 0.00054; %Wavelength in mm
Nyquest = lambda/(4*NA0); %Nyquest requirement for incoherent system for comparison
pixel = 0.003; %ideal pixel size in mm
mag = 10; %magnification of system
dx = pixel/mag; %Real space 'pixel' size in mm
N = 301; %desired # of pixels, make odd so we can see exact center
M = ceil(N/2);
fov = N*dx;
[xx,yy] = meshgrid([-N/2:N/2-1]*dx); %Spatial grid for use generating phase mask
DoF = 2*(lambda/(NA0^2)); %Native depth of field in this system

% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point

%% Generate Fresnel Disk
%Our ring is 0.45NA and is 1.3mm in diameter
%The fresnel ring equation is, for ring n, rn = sqrt((f+n*lambda/2)^2-f^2)
n = 0;
rn = 0;
fresnel = zeros(N,N);
while n<10 %Radius of miniscope
    n = n+1; %Index forward
    rn = sqrt(n*lambda*(f+(n*lambda)/2)); %Outer Radius of section
    %Convert radius to size in NA space
    fov = 1.3;%Size of aperture in mm
    NArn = rn/fov; %NA formula
    
    if mod(n,2) == 1
        fresnel(NAx.^2+NAy.^2 >= NArn.^2) = -pi;
    else
        fresnel(NAx.^2+NAy.^2 >= NArn.^2) = 0;
    end
end

imagesc(fresnel)