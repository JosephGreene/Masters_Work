%% Bessel Beam Analysis of EDoF Binary Axicon
%{
    Forth code designing an EDoF phase mask for incoherent/flourescent
    emission detection for neural tissue in mice. The goal of this code is
    to determine whether the beam produces is a bessel beam. This is
    determined by comparing axial and lateral resolution. A bessel beam
    elongates axial resolution without changing lateral. Else it is a
    spherical aberation.
%}
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

% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point
%% Loop to compare range of defocus of Axicon for different parameters
depth = 0.01;
step = 0.0001;
dz = -depth:step:depth; %small uniform field to check over, smaller dz
%3D Matrix that will hold 2D slices representing PSF at certain output
EDoF = zeros(N,N,length(dz));
FWHM = zeros(1,length(dz));

for n = 6.77
    %Generate mask and OTF/MTF
    if n == 0
        alpha = 0;
    else
        alpha =  NA0^2*f/(n*lambda)*10^-2;
    end
    binaxi = Generate_Binary_Axicon(alpha,NAx,NAy,0);
    %Generate on-axis MTF
    OTF_axi = acrr(binaxi);
    %Generate defocus
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF_axi = acrr(pre(binaxi.*defocus_prop));
        dMTF_axi = abs(dOTF_axi);
        iPSF = iF(dOTF_axi);
        EDoF(:,:,k) = iPSF;%iPSF = iF{OTF}
        
        %Find FWHM as function of dz
        maxValue = max(abs(iPSF(:)));
        % Find where it's more than half the max.
        aboveHalfMax = abs(EDoF(:,:,k)) > maxValue/2;
        % Get the first and last index where it's more than the half max.
        Index = find(aboveHalfMax);
        if length(Index) == 1
             FWHM(k) = 0.00005; %If just one pixel, FWHM is within peak pixel
        else
             [a,b] = ind2sub([301,301],Index); %Convert from linear indicies
             %Note, peak in center of bessel will ALWAYS be above FWHM
             FWHM(k)= max(sqrt(sum(([ceil(N/2),ceil(N/2)] - [a,b]).^2,2)))*0.0003; %Else find the maximum difference between last FWHM pixel and middle
        end
    end
    figure
    %Plot PSF for all dz @ y = 0 to determine depth of field created by this axicon
    %Center of Output is at N due to autocorrelation length being 2N-1
    imagesc(dz,M-50:M+50,squeeze(real(EDoF(M-50:M+50,M,:))))%squeeze compresses slice of 3d into 2d
    xlabel('Field depth, mm')
    ylabel('cross section of iPSF @ y = 0')
    title(['output @ y= 0|n = ' num2str(n)])
    colorbar
    
    figure
    plot(-depth:step:depth,FWHM);
    xlabel('Depth')
    ylabel('FWHM')
end
