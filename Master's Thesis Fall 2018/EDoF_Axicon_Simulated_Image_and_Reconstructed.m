%% Declare Setup and Aperature in NA space
%{
    Second code designing an EDoF phase mask for incoherent/flourescent
    emission detection for neural tissue in mice. The goal of this code is
    to determine the EDoF for different binary mask designs (of an axicon)
    in a noisey neural images. Also, later on, I look at different adjusted
    binaxicon masks that more approach fresnel disks. Also, I determine SBR 
    for each case.
%}
%Useful formulas
F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
acrr= @(x) iF(conj(F(x)).*F(x)); %autocorrelation using fourier transform (much faster)
pre = @(x) (x)/norm(x); %Preprocessing step to normalize input data for autocorrelation


%Physical parameters of miniscope (unchagable) in real space
NA0 = 0.45; %NA of the system -> range of ATF (2*NA is range of OTF)
f = 3.3; %Focus in mm
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

%Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point

%% Generate Simulated 3D Object of Defocues 'neurons'
%Generate binary axicon mask
n = 6.77; %Current param used in practice
alpha =  NA0^2*f/(n*lambda)*10^-2; %Slope of axicon
binaxi = Generate_Binary_Axicon(alpha,NAx,NAy,0);

%Generate 3D Object using waleed's code
dz = -0.5:0.0005:0.005; %Neural depth in mm
addpath('Waleeds_3D_Object')
[len, width] = size(dz);
simobj = sim_obj(N,N,width,8,[1]); %Only have one circle with radius 8*pixel/mag per dz slice

%Generate iPSF at every declared dz depth (adjusted for full neural depth)
%and defocued 'neruon' at each dz
%3D Matrix that will hold 2D slices representing iPSF at certain output
EDoF = zeros(N,N,length(dz));
convolved = zeros(size(EDoF));
%no_mask = zeros(size(EDoF));
for k = 1:length(dz)
    defocus = dz(k);%Defocus distance
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
    dOTF_axi = acrr(pre(binaxi.*defocus_prop)); %Defocus OTF
    EDoF(:,:,k) = iF(dOTF_axi);%iPSF = iF{OTF}, I only save EDoF for potential use
    convolved(:,:,k) = iF(dOTF_axi.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
end

%% Show Image and Initial Reconstruction
ref = sum((simobj),3); %Reference image, just each neuron compressed onto a 3D plane
figure
imagesc(ref)
title('Reference')

imager = sum((convolved),3); %Image received by full defocused field, noisy image
figure
imagesc(abs(imager))
title('imager reading')

on_axis_OTF = F(squeeze(EDoF(:,:,1001))); %Initial Reconstuction w/Tikhonov
rec = abs((iF((conj(on_axis_OTF).*F(imager))./(abs(on_axis_OTF).^2+0.006))));
figure
imagesc(rec)
title('reconstruction')

%% Optimal Reconstruction Parameter
u = 10.^[-5:0.25:5];

for i = 1:length(u)
    totale(i) = norm(((conj(on_axis_OTF).*F(imager))./(abs(on_axis_OTF).^2.+u(i)))-imager);
    approxe(i) = norm(((conj(on_axis_OTF).*ref)./(abs(on_axis_OTF).^2.+u(i)))-imager);
    noisee(i) = norm((conj(on_axis_OTF).*F(imager-ref))./(abs(on_axis_OTF).^2.+u(i)));
end

figure
hold on
plot(log10(u),totale)
plot(log10(u),approxe)
plot(log10(u),noisee)
hold off
xlabel('log10(u)')
ylabel('error')
title('Total Reference Object Reconstruction')
legend('Total error','Approx Error','Noise Error')

figure
plot(log10(u),approxe)
xlabel('log10(u)')
ylabel('Approx error')
title('Total Reference: Approx Error')
%% Now use a reference where only looking at neurons at depths +- 5 um (reconstruction range)
u = 10.^[-5:0.25:5];
newref = sum(simobj(:,:,991:end),3); %Sum last 20 images in simobj as reference
for i = 1:length(u)
    totale(i) = norm(((conj(on_axis_OTF).*F(imager))./(abs(on_axis_OTF).^2.+u(i)))-imager);
    approxe(i) = norm(((conj(on_axis_OTF).*ref)./(abs(on_axis_OTF).^2.+u(i)))-imager);
    noisee(i) = norm((conj(on_axis_OTF).*F(imager-ref))./(abs(on_axis_OTF).^2.+u(i)));
end

figure
hold on
plot(log10(u),totale)
plot(log10(u),approxe)
plot(log10(u),noisee)
hold off
xlabel('log10(u)')
ylabel('error')
title('Partial Reference Object Reconstruction')
legend('Total error','Approx Error','Noise Error')

figure
plot(log10(u),approxe)
xlabel('log10(u)')
ylabel('Approx error')
title('Total Reference: Approx Error')

%% Assume Background DC test
%Spectrum Analysis
figure
imagesc(abs(F(ref)))
title('Object Spectrum')
axis([110 190 110 190])

figure
imagesc(abs(F(imager)))
title('Noisy Object Spectrum')
axis([110 190 110 190])

%Remove DC
filt = ones(N,N);
filt(151,151) = 0;

figure
imagesc(abs(iF(F(ref).*filt)))
title('Ref no background')

filt = (abs(iF(F(imager).*filt)));
figure
imagesc(filt)
title('Noisy no background')

%% Attempt Background Removal Through Combined Hilbert Transform
%See Summer Report for details on combined hilbert transform
[n1,n2] = Signc2(N);
new_mask = binaxi.*n1.*n2; %Apply conditions for 2D combined hilbert transform
figure
imagesc(new_mask)
title('Hilbert Mask Design')

hilbert_convolved = zeros(size(EDoF));
for k = 1:length(dz)
    defocus = dz(k);%Defocus distance
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
    dOTF_axi = acrr(pre(new_mask.*defocus_prop)); %Defocus OTF
    hilbert_convolved(:,:,k) = iF(dOTF_axi.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
    if k == 1001
        new_on_axis = dOTF_axi; %At dz = 0, keep OTF
    end
end

%% Reconstruct New Mask
new_imager = sum((hilbert_convolved),3); %Image received by full defocused field, noisy image
figure
imagesc(abs(new_imager))
title('Hilbert mask imager reading')

u = 10.^[-6:0.1:6];
for i = 1:length(u)
    new_totale(i) = norm(((conj(new_on_axis).*F(new_imager))./(abs(new_on_axis).^2.+u(i)))-new_imager);
    new_approxe(i) = norm(((conj(new_on_axis).*ref)./(abs(new_on_axis).^2.+u(i)))-new_imager);
    new_noisee(i) = norm((conj(new_on_axis).*F(new_imager-ref))./(abs(new_on_axis).^2.+u(i)));
end

figure
hold on
plot(log10(u),new_totale)
plot(log10(u),new_approxe)
plot(log10(u),new_noisee)
hold off
xlabel('log10(u)')
ylabel('error')
title('New Mask Tikhonov Error')
legend('Total error','Approx Error','Noise Error')

figure
plot(log10(u),new_approxe)
xlabel('log10(u)')
ylabel('Approx error')
title('New Mask: Approx Error')

%% Visualize for several mu
figure
new_rec = abs((iF((conj(new_on_axis).*F(new_imager))./(abs(new_on_axis).^2+0.01))));
imagesc(new_rec)
title('Hilbert mask reconstruction: mu = 0.01')

figure
imagesc(abs((iF((conj(new_on_axis).*F(new_imager))./(abs(new_on_axis).^2+0.1)))))
title('Hilbert mask reconstruction: mu = 0.1')

figure
imagesc(abs((iF((conj(new_on_axis).*F(new_imager))./(abs(new_on_axis).^2+1)))))
title('Hilbert mask reconstruction: mu = 1')

%% Original Mask with HPF
HPF_mask = binaxi;
HPF_mask(M-10:M+10,M-10:M+10) = -1;
figure
imagesc(HPF_mask)
title('HPF Mask Design')


HPF_convolved = zeros(size(EDoF));
for k = 1:length(dz)
    defocus = dz(k);%Defocus distance
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
    dOTF_axi = acrr(pre(HPF_mask.*defocus_prop)); %Defocus OTF
    HPF_convolved(:,:,k) = iF(dOTF_axi.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
    if k == 1001
        HPF_on_axis = dOTF_axi; %At dz = 0, keep OTF
    end
end

%% Reconstruct HPF Mask
HPF_imager = sum((HPF_convolved),3); %Image received by full defocused field, noisy image
figure
imagesc(abs(HPF_imager))
title('HPF mask imager reading')

u = 10.^[-6:0.1:6];
for i = 1:length(u)
    HPF_totale(i) = norm(((conj(HPF_on_axis).*F(HPF_imager))./(abs(HPF_on_axis).^2.+u(i)))-HPF_imager);
    HPF_approxe(i) = norm(((conj(HPF_on_axis).*ref)./(abs(HPF_on_axis).^2.+u(i)))-HPF_imager);
    HPF_noisee(i) = norm((conj(HPF_on_axis).*F(HPF_imager-ref))./(abs(HPF_on_axis).^2.+u(i)));
end

figure
hold on
plot(log10(u),HPF_totale)
plot(log10(u),HPF_approxe)
plot(log10(u),HPF_noisee)
hold off
xlabel('log10(u)')
ylabel('error')
title('New Mask Tikhonov Error')
legend('Total error','Approx Error','Noise Error')

figure
plot(log10(u),new_approxe)
xlabel('log10(u)')
ylabel('Approx error')
title('New Mask: Approx Error')

%% Visualize for several mu
figure
HPF_rec = abs((iF((conj(HPF_on_axis).*F(HPF_imager))./(abs(HPF_on_axis).^2+0.01))));
imagesc(HPF_rec)
title('Hilbert mask reconstruction: mu = 0.01')

figure
imagesc(abs((iF((conj(HPF_on_axis).*F(HPF_imager))./(abs(HPF_on_axis).^2+0.1)))))
title('Hilbert mask reconstruction: mu = 0.1')

figure
imagesc(abs((iF((conj(HPF_on_axis).*F(HPF_imager))./(abs(HPF_on_axis).^2+1)))))
title('Hilbert mask reconstruction: mu = 1')
%% SBR For Original and New Mask
%Reference object for 'neurons' in +5um:-Num

for i = 0:0.5:50 %We start at dz = +5um
    subsim = simobj; 
    subsim(:,:,1:(length(dz)-i*2)) = 0;%Set undesired 'neurons' to zero
    subpixels = find(sum(subsim,3)); %find nonzeros indicies for 'imager' reading

    sbr_imager(i*2+1) = mean(abs(imager(subpixels)))./mean(abs(imager(:)));
    sbr_rec(i*2+1) = mean(rec(subpixels))./mean(rec(:));
    sbr_new_imager(i*2+1) = mean(abs(new_imager(subpixels)))./mean(abs(new_imager(:)));
    sbr_new_rec(i*2+1) = mean((new_rec(subpixels)))./mean(new_rec(:));
    sbr_filt(i*2+1) = mean((filt(subpixels)))./mean(filt(:));
    sbr_HPF_rec(i*2+1) = mean(HPF_rec(subpixels))./mean(HPF_rec(:));
end

%% Plot SBR
figure
hold on
plot(0:0.5:50,sbr_imager)
plot(0:0.5:50,sbr_rec)
plot(0:0.5:50,sbr_new_rec)
plot(0:0.5:50,sbr_filt)
plot(0:0.5:50,sbr_HPF_rec)
hold off
title('SBR Analysis')
ylabel('Average SBR')
xlabel('depth(um)-5(um)')
legend('Orig. Imag','Orig. Rec','Hilbert rec','Filtered', 'HPF rec')