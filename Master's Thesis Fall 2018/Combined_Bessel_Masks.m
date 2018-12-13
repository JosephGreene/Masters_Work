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
N = 301; %desired # of pixels, make odd so we can see exact center
M = ceil(N/2);
fov = N*dx; %N = 1000 in this case
[xx,yy] = meshgrid([-N/2:N/2-1]*dx); %Spatial grid for use generating phase mask
DoF = 2*(lambda/(NA0^2)); %Native depth of field in this system
%% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point

%% Combined Binary Plase Mask
n1 = 2;
n2 = 3;

alpha1 =  NA0^2*f/(n1*lambda)*10^-2;
alpha2 =  NA0^2*f/(n2*lambda)*10^-2;

binaxi1 = Generate_Binary_Axicon(alpha1,NAx,NAy,0);
binaxi2 = Generate_Binary_Axicon(alpha2,NAx,NAy,0);

figure
imagesc(binaxi1)
title(['binary axicon with n = ' num2str(n1)])
colorbar

figure
imagesc(binaxi2)
title(['binary axicon with n = ' num2str(n2)])
colorbar

combined = binaxi1.*binaxi2;
figure
imagesc(combined)
title('Combined Mask Design')
colorbar

%Generate defocus
dz = [-0.05:0.0005:0.05];
EDoF = zeros(N,N,length(dz));
for k = 1:length(dz)
    defocus = dz(k);%Defocus distance
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
    dOTF_axi = acrr(pre(combined.*defocus_prop));
    dMTF_axi = abs(dOTF_axi);
    EDoF(:,:,k) = iF(dOTF_axi);%iPSF = iF{OTF}
end

figure
%Plot PSF for all dz @ y = 0 to determine depth of field created by this axicon
%Center of Output is at N due to autocorrelation length being 2N-1
imagesc(dz,M-50:M+50,squeeze(real(EDoF(M-50:M+50,M,:))))%squeeze compresses slice of 3d into 2d
xlabel('Field depth, mm')
ylabel('cross section of iPSF @ y = 0')
title(['Combined Mask Output'])
colorbar