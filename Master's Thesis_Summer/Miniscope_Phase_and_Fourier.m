%% Mathematically Simulating a 4F System
%{
In this following report I will utilize the phase/propegation model and
fourier model to simulate how light travels through a 4F system. I develop
this code mostly as a reference guide for future-more complex-reports.

For the purpose of this code, I will assume uniform focal distances for
each lens, the mask is between the two lenses in the fourier plane and 
%}

%% Phase description -> Assumes Coherence
%Define Parameters
N = 400; %N by N images
M = floor(N/2);
delta = 2.2; %Pixel Pitch in microns
ds = 10; %Downsampling factor (for generating waveform)
lambda = 0.5; %Wavelength light
f = 3300; %Focal Length um
fov = 1/(N*delta);
[xx,yy] = meshgrid([-M:M-1]*delta);
[uu,vv] = meshgrid([-M:M-1]*fov);
k = 2*pi/lambda;
k2 = pi*lambda*(uu.^2+vv.^2);

I1 = load('I1.mat');
orig = I1.I1;

%Plot Image and used object
figure;
subplot(1,2,1)
imagesc(orig);
title('original image')
axis('square')

o1 = imresize(orig,[M,M]);
o1 = padarray(o1,[M/2,M/2],'both');

subplot(1,2,2)
imagesc(o1)
title('padded image')
axis('square')

%Progegation Equations
F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
prop = @(f,h) iF(F(f).*h);
df = exp(1j*k2*f);

%Lens function
lens = exp(1i*k/2/f.*(xx.^2+yy.^2));  %Lens Transmittance Func
input = zeros(N,N);
input(M+1,M+1) = 1; %Input Delta

field_lens1 = prop(o1,df).*lens; %Field after first lens
%field_lens1 = ones(N,N); %Ideal Field with delta input
mask1 = wrapTo2Pi(2*pi/lambda*(1.5-1)*6.5*10^7*(xx.^3+yy.^3));%Phase function of cubic phase mask
mask1((xx.^2+yy.^2) >= 0.95*N^2) = 0;% Circular apature function 
mask2 = ones(N,N); %Ideal Apature
field_mask = prop(field_lens1,df).*mask2; %Field after mask
field_lens2 = prop(field_mask,df).*lens; %Field after second lens
output = prop(field_lens2,df); %Output

figure
subplot(1,2,1)
imagesc(abs(field_mask))
title('Fourier Plane: Post Mask')
axis([M-30 M+30 M-30 M+30])
axis('square')
subplot(1,2,2)
imagesc(abs(output))
title('Output')
axis('square')

%% Fourier Optics -> Incoherent Source
%{
For the described 4f system:

PSF = F(mask)
ATF = iF(F(mask)) = mask
OTF = autocorrelation(ATF,ATF)
Output = iF(F(abs(o(x,y)).^2).*OTF)

%%Note Images loaded from the raw data represents abs(o(x,y)).^2, i.e
intensity measurements
%}

%Load raw data
raw = VideoReader('raw1.mov');

%First Practice Spot (Away from central/high noise)
raw.CurrentTime = 10.6; %Spot Upper Left
x = [round(0.55*N):round(0.57*N)]; %x coodinates
y = [round(0.1375*N):round(0.17*N)]; %y coodinates

%Second Practice Image (In high noise in center)
% raw.CurrentTime = 6.9; %Center Spot
% x = [round(0.455*N):round(0.505*N)]; %x coordinates
% y = [round(0.45*N):round(0.5175*N)]; %y coodinates

frame = imresize(im2double(rgb2gray(readFrame(raw))),[N,N]); 
% frame = padarray(frame,[M/2,M/2],'both');

%% Use of Phase Objects for HP filtering -> Coherent
%{
This section discusses the implementation of a spiral phase mask (spm) as
described by "Spiral Phase Contrast Imaging in Microscopy". Note this
technique performs high pass filtering of an input object, but only works
on phase objects. As an experiment, I will simulate the effects of a phase
mask on an object utilized by the paper as well as an image from the
miniscope. The miniscope is an amplitude image - so this will never work,
but I'd like to observe the 'ideal case' for a HP filtering mask.
%}
o2 = imread('5_obj.png');
o2 = imresize(im2double(rgb2gray(o2)),[N,N]);
spm = Spiral_Phase_Mask(1,N);
output1 = iF(F(exp(1j*frame)).*spm);
output2 = iF(F(exp(1j*o2)).*spm);

figure
subplot(1,2,1)
imagesc(frame)
title('Orig')
axis('square')
subplot(1,2,2)
imagesc(abs(output1))
title('SPM Filtered')
axis('square')

figure
subplot(1,2,1)
imagesc(o2)
title('Orig')
axis('square')
subplot(1,2,2)
imagesc(abs(output2))
title('SPM Filtered')
axis('square')

%% Filtering a video frame by frame
%{
Filtere each frame of original movie by spm and create output movie
%}

raw.CurrentTime = 0;
filt = VideoWriter('filtered_new.avi');
open(filt)
while(readFrame(raw))
    new_frame = iF(F(exp(1j*imresize(im2double(rgb2gray(readFrame(raw))),[N,N]))).*spm);
    writeVideo(filt,new_frame./max(new_frame(:)))
end
close(filt)

%% 
