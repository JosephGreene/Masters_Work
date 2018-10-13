%% Binary Masks EDoF properties
F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)
pre = @(x) (x)/norm(x); %Preprocessing step to normalize input data for autocorrelation


%Physical parameters of miniscope (unchagable) in real space
NA0 = 0.45; %NA of the system -> range of ATF (2*NA is range of OTF)
f = 3.3; %mm
dimag = 10; %mm
lambda = 0.00054; %Wavelength in mm
Nyquest = lambda/(4*NA0); %Nyquest requirement for incoherent system for comparison
pixel = 0.003; %ideal pixel size in mm
mag = 10; %magnification of system
dx = pixel/mag; %Real space 'pixel' size in mm
N = 401; %desired # of pixels, make odd so we can see exact center
fov = N*dx; %N = 1000 in this case
[xx,yy] = meshgrid([-N/2:N/2-1]*dx); %Spatial grid for use generating phase mask
DoF = 2*(lambda/(NA0^2)); %Native depth of field in this system
%% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point

%% Define Nonspherical Masks

%[1] tangent mask
a = 37.53*10^-2;
b = 1.27;
%[2] Binarized Cubic Phase Mask
a2 = 6.5*10^2;

tangent = wrapTo2Pi(a*NAx.^2*tan(b*NAx)+a*yy.^2*tan(b*NAy));
tangent((NAx.^2+NAy.^2) >= (0.45)^2) = 0; %Mimic Aperture
%Visualize masks
figure
imagesc(tangent)
title('Tangent Mask')

cubic = a2*(NAx.^3+NAy.^3);
cubic((NAx.^2+NAy.^2) >= (0.45)^2) = 0; %Mimic Aperture
wrapped = wrapTo2Pi(cubic);
bincubic = wrapped;
bincubic(bincubic >= pi) = 2*pi;
bincubic(bincubic < pi) = 0;
figure
imagesc(bincubic)
title('Binary Cubic Mask')
%% iPSFs for defined masks
dz = [-0.015:0.00025:0.015];
%3D Matrix that will hold 2D slices representing PSF at certain output
EDoF_tan = zeros(N,N,length(dz));
EDoF_cub= zeros(N,N,length(dz));

%Generate OTF

OTF_tan = acrr(tangent);
OTF_cub = acrr(bincubic);

%Generate iPSF with defocus
for k = 1:length(dz)
    defocus = dz(k);%Defocus distance
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
    dOTF_tan = acrr(pre(tangent.*defocus_prop));
    dOTF_cub = acrr(pre(bincubic.*defocus_prop));
    EDoF_tan(:,:,k) = iF(dOTF_tan);%iPSF = iF{OTF}
    EDoF_cub(:,:,k) = iF(dOTF_cub);%iPSF = iF{OTF}
end
figure
%Plot PSF for all dz @ y = 0 to determine depth of field created by this axicon
%Center of Output is at N due to autocorrelation length being 2N-1
imagesc(dz,150:250,squeeze(real(EDoF_tan(150:250,floor(N/2)+1,:))))%squeeze compresses slice of 3d into 2d
xlabel('Field depth, mm')
ylabel('cross section of iPSF @ y = 0')
title(['output @ y= 0 for tan'])
colorbar

figure
imagesc(dz,150:250,squeeze(real(EDoF_cub(150:250,floor(N/2)+1,:))))%squeeze compresses slice of 3d into 2d
xlabel('Field depth, mm')
ylabel('cross section of iPSF @ y = 0')
title(['output @ y= 0 for binary cubic'])
colorbar
