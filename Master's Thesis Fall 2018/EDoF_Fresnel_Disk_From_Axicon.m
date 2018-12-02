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
[xx,yy] = meshgrid([-M:M-1]*dx); %Spatial grid for use generating phase mask
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

%%
%Generate 3D Object using waleed's code
dz = -0.1:0.0005:0.005; %Neural depth in mm -> note I reduce total depth to speed up code form prev ex
addpath('Waleeds_3D_Object')
[len, width] = size(dz);
simobj = sim_obj(N,N,width,2,[2]); %Only have one circle with radius 8*pixel/mag per dz slice

%% Fresnel Disk 1: Binary Axicon with Inner circle set to pi shift
%Generate iPSF at every declared dz depth (adjusted for full neural depth)
%and defocued 'neruon' at each dz
%3D Matrix that will hold 2D slices representing iPSF at certain output
no_mask = zeros(N,N);
noMask = zeros(N,N);
no_mask((NAx.^2+NAy.^2) <= 0.45.^2) = 1;

radii = [0,0.05,0.075,0.314]; %Radius of inner ring increasing at geometric average

for j = 1:length(radii) %Adjust radius in loop
    image = zeros(N,N);
    new_design = binaxi;
    new_design(NAx.^2 + NAy.^2 <= radii(j).^2) = -1;
    figure
    imagesc(new_design)
    title(['New Disk w/ inner radius of NA = ' num2str(radii(j))])
   
    %Determine output image
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF = acrr(pre(new_design.*defocus_prop)); %Defocus OTF
        image = image + iF(dOTF.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        if dz(k) == 0
            on_axis = dOTF;
        end
        %For only 1 iteration, save an image of just the defocused neural
        %tissue (no mask, baseline)
        if j == 1
            noMask = noMask + iF(acrr(pre(no_mask.*defocus_prop)).*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        end
    end

    %Reconstruction
    rec = abs((iF((conj(on_axis).*F(image))./(abs(on_axis).^2+0.01))));
    figure
    imagesc(rec)
    title(['Fresnel Disk Rec w/ inner radius of NA = ' num2str(radii(j))])
    
    %Determine SBR for each case
    for i = 0:0.5:25 %We start at dz = +5um
        subslice = simobj(:,:,(length(dz)-i*2)); %get slice at depth
        subpixels = find(subslice); %find nonzeros indicies for 'imager' reading
        invpixels = find(~subslice);
        
        sbr_image(i*2+1,j) = mean(abs(image(subpixels)))./mean(abs(image(invpixels)));
        sbr_rec(i*2+1,j) = mean(abs(rec(subpixels)))./mean(abs(rec(:)));
        if j == 1
           sbr_mask(i*2+1) = mean(abs(noMask(subpixels)))./mean(abs(noMask(invpixels)));
        end
    end
end

figure
imagesc(no_mask)
title('Aperture')

figure
imagesc(abs(noMask))
title('Output Image: No Mask')
    
%% Plot SBR for important cases
figure('Name','SBR Fresnel Disk')
hold on
plot(5:-0.5:-20,sbr_mask)
plot(5:-0.5:-20,sbr_image(:,1))
legendcell{1} = 'No Mask';
legendcell{2} = 'Original Image';
for i = 1:j
    plot(5:-0.5:-20,sbr_rec(:,i))
    legendcell{i+2} = (['rec NA = ' num2str(radii(i))]);
end
title('SBR: Fresnel Disk')
xlabel('Depth um')
ylabel('SBR')
set(gca,'Xdir','reverse')
legend(legendcell)

%% Mask Based on Learned phase coded aperture for the benefit of depth of field extension using phase
phase = [0.2,0.43,0.5,0.6,1]*2*pi; %Phase Shift
r_ratio = 0.8; %Ratio of outer and inner ratio

%Test Dependence on ratio
for j = 1:length(phase)
    coded = zeros(N,N);
    coded((NAx.^2+NAy.^2) <= (r_ratio*NA0).^2) = exp(1i*0);
    coded((NAx.^2+NAy.^2) > (r_ratio*NA0).^2) = real(exp(1i*phase(j)));
    coded((NAx.^2+NAy.^2) > (NA0).^2) = 0;
    
    figure
    imagesc(coded)
    title(['Coded Mask w/ Phase = ' num2str(phase(j))])
    
    %Determine output image
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF = acrr(pre(coded.*defocus_prop)); %Defocus OTF
        image = image + iF(dOTF.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        if dz(k) == 0
            on_axis = dOTF;
        end
    end
    
    %Reconstruction
    rec = abs((iF((conj(on_axis).*F(image))./(abs(on_axis).^2+0.1))));
    figure
    imagesc(rec)
    title(['Coded Mask Rec w/ Phase = ' num2str(phase(j))])
    
    %Determine SBR for each case
    for i = 0:0.5:25 %We start at dz = +5um
        subslice = simobj(:,:,(length(dz)-i*2)); %get slice at depth
        subpixels = find(subslice); %find nonzeros indicies for 'imager' reading
        invpixels = find(~subslice);
        
        sbr_image(i*2+1,j) = mean(abs(image(subpixels)))./mean(abs(image(invpixels)));
        sbr_rec(i*2+1,j) = mean(abs(rec(subpixels)))./mean(abs(rec(:)));
    end
end
%% Plot SBR for important cases
figure('Name','SBR Height')
hold on
plot(5:-0.5:-20,sbr_mask)
legendcell{1} = 'No Mask';
for i = 1:j
    plot(5:-0.5:-20,sbr_rec(:,i))
    legendcell{i+1} = (['rec phase = ' num2str(phase(i))]);
end
title('SBR: Adjust Phase')
xlabel('Depth um')
ylabel('SBR')
set(gca,'Xdir','reverse')
legend(legendcell)

%% depth of field extension adjusting ratio
phase = pi; %Phase Shift
r_ratio = 0.6:0.1:0.8; %Ratio of outer and inner ratio

%Test Dependence on ratio
for j = 1:length(r_ratio)
    coded = zeros(N,N);
    coded((NAx.^2+NAy.^2) <= (r_ratio(j)*NA0).^2) = exp(1i*0);
    coded((NAx.^2+NAy.^2) > (r_ratio(j)*NA0).^2) = real(exp(1i*phase));
    coded((NAx.^2+NAy.^2) > (NA0).^2) = 0;
    figure
    imagesc(coded)
    title(['Coded Mask w/ Radius Ratio = ' num2str(r_ratio(j))])
    
    %Determine output image
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF = acrr(pre(coded.*defocus_prop)); %Defocus OTF
        image = image + iF(dOTF.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        if dz(k) == 0
            on_axis = dOTF;
        end
    end
    
    %Reconstruction
    rec = abs((iF((conj(on_axis).*F(image))./(abs(on_axis).^2+0.1))));
    figure
    imagesc(rec)
    title(['Coded Mask Rec w/ Radius Ratio = ' num2str(r_ratio(j))])
    
    %Determine SBR for each case
    for i = 0:0.5:25 %We start at dz = +5um
        subslice = simobj(:,:,(length(dz)-i*2)); %get slice at depth
        subpixels = find(subslice); %find nonzeros indicies for 'imager' reading
        invpixels = find(~subslice);
        
        sbr_image(i*2+1,j) = mean(abs(image(subpixels)))./mean(abs(image(invpixels)));
        sbr_rec(i*2+1,j) = mean(abs(rec(subpixels)))./mean(abs(rec(:)));
    end
end
%% Plot SBR for important cases
figure('Name','SBR r ratio')
hold on
plot(5:-0.5:-20,sbr_mask)
legendcell{1} = 'No Mask';
for i = 1:j
    plot(5:-0.5:-20,sbr_rec(:,i))
    legendcell{i+1} = (['rec r ratio = ' num2str(r_ratio(i))]);
end
title('SBR: Adjust Radius')
xlabel('Depth um')
ylabel('SBR')
set(gca,'Xdir','reverse')
legend(legendcell)

%% Determine effectiveness of a design with two concentric rings
phase = pi; %Phase Shift
dNA = NAx(1,2) - NAx(1,1); %Stepsize in NA space
r1 = 0.45; %Radius first ring
t1 = 7; %Thickness first ring
r2 = 0.4;
t2 = 7:3:16;

%Test Dependence on inner ring radius
for j = 1:length(t2)
    coded = zeros(N,N);
    %Outer Ring,
    coded((NAx.^2+NAy.^2) <= (r1).^2) = exp(1i*0);
    coded((NAx.^2+NAy.^2) < (r1-t1*dNA).^2) = real(exp(1i*phase));

    %Inner Ring

    coded((NAx.^2+NAy.^2) <= (r2).^2) = exp(1i*0);
    coded((NAx.^2+NAy.^2) < (r2-t2(j)*dNA).^2) = real(exp(1i*phase));

    coded((NAx.^2+NAy.^2) > (NA0).^2) = real(exp(1i*phase));

    figure
    imagesc(coded)
    title(['Coded Mask w/ Inner Radius = ' num2str(t2(j))])
    
    %Determine output image
    for k = 1:length(dz)
        defocus = dz(k);%Defocus distance
        defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2)); %Fresnel Kernel
        dOTF = acrr(pre(coded.*defocus_prop)); %Defocus OTF
        image = image + iF(dOTF.*F(squeeze(simobj(:,:,k))));%Squeeze converts 3D matrix slice into 2D matrix 
        if dz(k) == 0
            on_axis = dOTF;
        end
    end
    
    %Reconstruction
    rec = abs((iF((conj(on_axis).*F(image))./(abs(on_axis).^2+0.1))));
    figure
    imagesc(rec)
    title(['Coded Mask Rec w/ Inner Radius = ' num2str(t2(j))])
    
    %Determine SBR for each case
    for i = 0:0.5:25 %We start at dz = +5um
        subslice = simobj(:,:,(length(dz)-i*2)); %get slice at depth
        subpixels = find(subslice); %find nonzeros indicies for 'imager' reading
        invpixels = find(~subslice);
        
        sbr_image(i*2+1,j) = mean(abs(image(subpixels)))./mean(abs(image(invpixels)));
        sbr_rec(i*2+1,j) = mean(abs(rec(subpixels)))./mean(abs(rec(:)));
    end
end

%% Plot SBR for important cases
figure('Name','SBR r ratio')
hold on
plot(5:-0.5:-20,sbr_mask)
legendcell{1} = 'No Mask';
for i = 1:j
    plot(5:-0.5:-20,sbr_rec(:,i))
    legendcell{i+1} = (['rec t2 = ' num2str(t2(i))]);
end
title('SBR: Adjust Radius')
xlabel('Depth um')
ylabel('SBR')
set(gca,'Xdir','reverse')
legend(legendcell)