%% Binary Axicon EDoF properties
%{
    First code designing an EDoF phase mask for incoherent/flourescent
    emission detection for neural tissue in mice. The goal of this code is
    to determine the EDoF for different binary mask designs (of an axicon)
    in a noiseless chunk of neural tissue.
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
%% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point
%% Loop to compare range of defocus of Axicon for different parameters
dz = [-0.05:0.0005:0.05];
ip_axi = zeros(1,length(dz));
%3D Matrix that will hold 2D slices representing PSF at certain output
EDoF = zeros(N,N,length(dz));

for n = 4
    %Generate mask and OTF/MTF
    if n == 0
        alpha = 0;
    else
        alpha =  NA0^2*f/(n*lambda)*10^-2;
    end
    binaxi = Generate_Binary_Axicon(alpha,NAx,NAy,0);
    %Generate on-axis MTF
    OTF_axi = acrr(binaxi);
    MTF_axi = abs(OTF_axi);
    figure
    imagesc(binaxi)
    title(['binary axicon with n = ' num2str(n)])
    colorbar
    figure
    imagesc(OTF_axi);
    title(['binary axicon OTF with n = ' num2str(n)])
    colorbar
    %Generate defocus
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF_axi = acrr(pre(binaxi.*defocus_prop));
        dMTF_axi = abs(dOTF_axi);
        EDoF(:,:,k) = iF(dOTF_axi);%iPSF = iF{OTF}
    end
    figure
    %Plot PSF for all dz @ y = 0 to determine depth of field created by this axicon
    %Center of Output is at N due to autocorrelation length being 2N-1
    imagesc(dz,M-50:M+50,squeeze(real(EDoF(M-50:M+50,M,:))))%squeeze compresses slice of 3d into 2d
    xlabel('Field depth, mm')
    ylabel('cross section of iPSF @ y = 0')
    title(['output @ y= 0|n = ' num2str(n)])
    colorbar

    %plot(ip_axi./max(ip_axi))
end

%% Generate 3D object to use with iPSF generated -> Waleed's code
addpath('Waleeds_3D_Object')
[len, width] = size(dz);
simobj = sim_obj(N,N,width,8,[1]);

%% Generate Output!
%Recall imager is 2D so the imager will see the convolution in x,y and
%projection in z. I.E summation of 2D convolution
convolved = zeros(size(EDoF));
[x,y,z] = size(convolved);
for i = 1:z
    convolved(:,:,i) = iF(F(squeeze(EDoF(:,:,i))).*F(squeeze(simobj(:,:,i)))); 
end

imager = sum((convolved),3);
figure
imagesc(abs(imager))
title('imager reading')

% ref = sum(simobj,3);
% figure
% imagesc(ref)
% title('reference')
% 
on_axis_OTF = F(squeeze(EDoF(:,:,1001)));
figure
imagesc(abs((iF((conj(on_axis_OTF).*F(imager))./(abs(on_axis_OTF).^2+0.006)))))
title('reconstruction')
%% Approximation Error
approxe = ones(1,7);
for u = 10.^[-3:3]
    i = 1;
    approxe(i) = norm(conj(on_axis_OTF).*F(imager)./(abs(on_axis_OTF).^2+u).*ref-ref);
    i= i+1;
end

%% Bonus Round Combining Masks
dz = [-0.02:0.0005:0.02];
%3D Matrix that will hold 2D slices representing PSF at certain output
EDoF2 = zeros(N,N,length(dz));

for n2 = 2
    %Generate mask and OTF/MTF
    if n2 == 0
        alpha2 = 0;
    else
        alpha2 =  NA0^2*f/(n2*lambda)*10^-2;
    end
    binaxi2 = Generate_Binary_Axicon(alpha2,NAx,NAy,0);
    %Generate on-axis MTF
    OTF_axi2 = acrr(binaxi2);

    figure
    imagesc(binaxi2.*binaxi)
    title(['binary axicon with n = ' num2str(n)])
    colorbar
    figure
    imagesc(OTF_axi2);
    title(['binary axicon OTF with n = ' num2str(n)])
    colorbar
    %Generate defocus
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF_axi2 = acrr(pre(binaxi2.*binaxi.*defocus_prop));
        EDoF2(:,:,k) = iF(dOTF_axi2);%iPSF = iF{OTF}
    end
    figure
    %Plot PSF for all dz @ y = 0 to determine depth of field created by this axicon
    %Center of Output is at N due to autocorrelation length being 2N-1
    imagesc(dz,M-50:M+50,squeeze(real(EDoF2(M-50:M+50,M,:))))%squeeze compresses slice of 3d into 2d
    xlabel('Field depth, mm')
    ylabel('cross section of iPSF @ y = 0')
    title(['output @ y= 0|n = ' num2str(n)])
    colorbar

    %plot(ip_axi./max(ip_axi))
end