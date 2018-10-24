%% Declare Setup and Aperature in NA space
%{
    Third code designing an EDoF phase mask for incoherent/flourescent
    emission detection for neural tissue in mice. The goal of this code is
    to determine the a qualitatively satisfying fresnel disk based on a
    previous binary axicon design that best 'kills' the ~equipower rings
    created by bessel functions while retaining a decent SBR
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

%% Generate orignal mask and 3D object
%Generate binary axicon mask
n = 6.77; %Current param used in practice
alpha =  NA0^2*f/(n*lambda)*10^-2; %Slope of axicon
binaxi = Generate_Binary_Axicon(alpha,NAx,NAy,0);
figure
imagesc(binaxi)
title('original mask')

%Generate 3D Object using waleed's code
dz = -0.5:0.0005:0.005; %Neural depth in mm -> note I reduce total depth to speed up code form prev ex
addpath('Waleeds_3D_Object')
[len, width] = size(dz);
simobj = sim_obj(N,N,width,8,[1]); %Only have one circle with radius 8*pixel/mag per dz slice

%% Fresnel Disk 1: Binary Axicon with Inner circle set to pi shift
%Generate iPSF at every declared dz depth (adjusted for full neural depth)
%and defocued 'neruon' at each dz
%3D Matrix that will hold 2D slices representing iPSF at certain output
no_mask = zeros(N,N);
radii = [0,0.0314,0.05,0.1,0.314]; %Radius of inner ring increasing at geometric average

for j = 1:length(radii) %Adjust radius in loop
    image = zeros(N,N);
    fresnel = binaxi;
    fresnel(NAx.^2 + NAy.^2 <= radii(j).^2) = -1;
    figure
    imagesc(fresnel)
    title(['Fresnel Disk w/ inner radius of NA = ' num2str(radii(j))])
   
    %Determine output image
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF = acrr(pre(fresnel.*defocus_prop)); %Defocus OTF
        image = image + iF(dOTF.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        if dz(k) == 0
            on_axis = dOTF;
        end
        %For only 1 iteration, save an image of just the defocused neural
        %tissue (no mask, baseline)
        if j == 1
            no_mask = no_mask + iF(acrr(defocus_prop).*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        end
    end
    figure
    imagesc(abs(image))
    title('Output Image')
    
    %Reconstruction
    rec = abs((iF((conj(on_axis).*F(image))./(abs(on_axis).^2+0.01))));
    figure
    imagesc(rec)
    title(['Fresnel Disk Rec w/ inner radius of NA = ' num2str(radii(j))])
    
    %Determine SBR for each case
    for i = 0.5:0.5:50 %We start at dz = +5um
        subsim = simobj; 
        subsim(:,:,1:(length(dz)-i*2)) = 0;%Set undesired 'neurons' to zero
        subpixels = find(sum(subsim,3)); %find nonzeros indicies for 'imager' reading

        sbr_image(i*2+1,j) = mean(abs(image(subpixels)))./mean(abs(image(:)));
        sbr_rec(i*2+1,j) = mean(rec(subpixels))./mean(rec(:));
        if j == 1
           sbr_mask(i*2+1) = mean(abs(no_mask(subpixels)))./mean(abs(no_mask(:)));
        end
    end
end
figure
imagesc(abs(no_mask))
title('Maskless Image')

%% Plot SBR for important cases
figure
hold on
plot(sbr_mask)
plot(sbr_image(:,1))
legendcell{1} = 'No Mask';
legendcell{2} = 'Original Image';
for i = 1:j
    plot(sbr_rec(:,i))
    legendcell{i+2} = (['rec NA = ' num2str(radii(i))]);
end
title('SBR')
legend(legendcell)