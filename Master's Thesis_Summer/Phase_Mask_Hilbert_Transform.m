%% Phase Mask Design: Optical Hilbert Transform
%{
In this following report, I will use MATLAB to improve images of neural
activity in small animals to design a better method of collecting data in
the Boston University Freedom Scope. 

Below, I will define several parameters, common venacular, and abreviations
I will refer to in this project:
    *F{} = Fourier Transform

    *iF{} = Inverse Fourier Transform

    *Coherent Light: Light that exhibits consistent phase relationships
    between electric field values at different locations (Spatial
    Coherence) or times (Temporal Coherence). Multiple coherent light
    sources have the same frequency and may produce interference pattens.

    *Incoherent Light: Light produced by emitted photons of different
    frequencies. Unlike Coherent sources, interference patterns are
    indistinquishable making it impossible to measure the electric field.
    Instead incoherent imaging systems must measure differences in
    intensity (E^2, so randomly generated, muddied phase term eliminated).

    *Focal Plane: The location in front and in back of an imaging system
    where light is optimally focused when exiting. A point source
    generating light on the front focal plane appears in-focus when imaged
    at the back focal plane.

    *cPSF: coherent Point Spread Function, or the response of an imaging
    system when imaging a coherent point source. As a note, this describes
    the response in the electric field generated by the point source. An
    image may be approximated as a superposition of point sources so the
    output of an optical system may be calculated by the 2D convolution
    between the input image and cPSF

    *iPSF: incoherent Point Spread Function, or the response of an imaging
    system when imaging a incoherent point source. As a note, this describes
    the response in the intensity generated by the point source. Since
    intensity = electric_field^2, iPSF = |cPSF|^2.

    *ATF: Amplitude Transfer Function, recalling that a lens performs a
    fourier transform of in focus light, the amplitude transfer function
    represents how the point spread function appears in frequency, or the
    fourier transform of the cPSF. Computationally, I may also compute the
    output of an optical system as iF{F{o1}.*ATF}, for some input object
    o1. This is computationally more efficient than calculating the
    convolution because matlab's fourier transform is ~nln(n) mutilplications
    for an nxn object, while the convolution is ~n^2.

    *OTF: Optical Transfer Function, or the fourier domain representation
    of the iPSF. Since the ATF is defined as F{cPSF}, the OTF is F{iPSF =
    |cPSF|^2} = a(ATF,ATF), where A represent the autocorrelation.

    *4f system: An optical system that comprises of 2 lenses. In practice,
    the space between the two lenses is the fourier domain spectrum of the
    input light, allowing researchers to perform filtering by placing masks
    in between the lenses. 

    *Phase Mask: A mask that adjust the phase of incoming light to generate
    certain interference spectrum on the imager. Computationally, this
    allows us to remove noise, or blur incoming light by a uniform PSF,
    allowing us to computationally deblur the output image for a sharper
    focus.

    *Miniscope: a minaturized microscope designed to fasten on the heads of
    small animals and transmit images of neural activity.

The goal of this project is to produce a phase mask which will remove
background noise present in images collected by Boston University's
miniscope, the freedom scope. The freedom scope is a miniaturized 4f system
which monitors a incoherent flourescent staining agent designed to emit
flourophores whenever a stained neuron fires. This allows the researcher to
individual neural activity as the animals interact with the world. I plan
on doing this by researching phase mask designs to see which ones could
produce an OTF that performs a high pass filtering on the generated image
to remove the nonvarying background term in reconstruction.
%}

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
%% Generalized Cubic Phase Mask
%{
This phase mask design is popular in extended depth of field applications
(blurring points on and near the focal plane in the same way so we may
computationally deblur the image to have a larger portion in focus) [1]. I
will investigate how it's OTF compares to our desired parameter when it is
a phase wrapped, continuous function as well as a binary representation.
%}

%Define Parameters
n_mat = 1.5; %Index of refraction for material
lambda = 500*10^-9; %normalized to m, blue light
k = 2*pi/lambda;
p = 2.5*10^-6; %pixel pitch in m of miniscope
f = 3.3*10^-3; %Focal length
DS = 10; %Downsample factor (How much less sampling is then pixel pitch)
delta = p/DS; %sampling rate in spatial domain (to prevent aliasing in generated field) 
N = 400; %Number of spatial domain samples
M = floor(N/2);
Mask = 0.95*M*delta; %Fraction of mask unobscured by apature

[xx,yy] = meshgrid([-M:M-1]*delta); %Spatial Sampling
[uu,vv] = meshgrid([-M:M-1]*1/(2*N*delta)); %S. Frequency Sampling

%Generate Mask
gcpm = wrapTo2Pi(2*pi/lambda*(n_mat-1)*6.5*10^7*(xx.^3+yy.^3));%Phase function of cubic phase mask
gcpm((xx.^2+yy.^2) >= Mask^2) = 0;% Circular apature function

GCPM = F(gcpm);

OTF = xcorr2(GCPM,GCPM);
figure;
imagesc(abs(gcpm))
title('Phase Mask Design')

figure;
imagesc(abs(GCPM))
title('ATF')

figure;
imagesc(abs(OTF))
title('OTF')

%% Inverse Rect wave
rect = zeros(1,1000);
rect(490:510) = 1;
figure;
plot(abs(F(rect)))

irect = ones(1,1000);
irect(490:510) = 0;
figure;
plot(abs(F(irect)))

%% Attempt at Inverse Rect Phase mask
q = 2;
d1 = 5;
d2 = 20;
irec_mask = zeros(N,N);
irec_mask(M-d2:M+d2,M-d2:M+d2) = 2*pi/lambda*(n_mat-1);
irec_mask(M-d1:M+d1,M-d1:M+d1) = q*pi/lambda*(n_mat-1);
figure;
imagesc(irec_mask)
figure
imagesc(abs(F(irec_mask)))
colorbar

%% Spiral Phase Mask (Continuous)
%{
In this section, I analyze the fourier domain characteristics of a spiral
phase mask proposed in "Spiral phase mask shadow-imaging for 3D measurement
of flow fields". To generate the phase mask, I treat each ring as wrapped
exponential decay functions afterwords normalized by the desired phase
delay value. While there are other ways to accomplish this task, I found
this method to be the easiest to execute.
%}

q = 4; %Number of divisions
N = 400; %Number of spatial domain samples
M = floor(N/2);

[xx,yy] = meshgrid([-M:M-1]); %Spatial Sampling
spiral = zeros(size(xx)); %Make Zero Array
spiral((xx.^2+yy.^2) < 200.^2) = 1;
for n = 1:q
    if(n==1) %Inner spiral is twice as large as secondary ones
        in = 0;
        out = round(M*2/(q+1));
        theta = wrapTo2Pi(atan2(-yy,-xx));
    else %Generate secondary spirals
        in = n*round(M/(q+1));
        out = (n+1)*round(M/(q+1));
        theta = wrapTo2Pi((2*n-1)*atan2(-yy,-xx));
    end
    spiral((xx.^2+yy.^2)<=out.^2 & (xx.^2+yy.^2)>=in.^2) = exp(-theta((xx.^2+yy.^2)<=out.^2 & (xx.^2+yy.^2)>=in.^2));
end
spiral = spiral./max(spiral(:));
figure
imagesc(real(spiral));
colorbar

figure;
imagesc(abs(F(spiral)))

%% Spiral Phase Mask (Imaginary exponential)
%{
This phase mask design utilizes the structure outlined in the previous
section, however, I ensure the decay is symmetry respective to both sides
(I.E, the exponential is imaginary). As a result, while the previous mask
did not filter DC noise, this one does. This result is intuitive because an
imaginary exponential is decomposed into pure sinusoids, so no DC
component should arrise in the fourier plane.
%}
q = 1; %Number of divisions
N = 400; %Number of spatial domain samples
M = floor(N/2);

[xx,yy] = meshgrid([-M:M-1]); %Spatial Sampling
spiral = zeros(size(xx)); %Make Zero Array
spiral((xx.^2+yy.^2) < 200.^2) = 1;
for n = 1:q
    if(n==1) %Inner spiral is twice as large as secondary ones
        in = 0;
        out = round(M*2/(q+1));
        theta = wrapTo2Pi(atan2(-yy,-xx));
    else %Generate secondary spirals
        in = n*round(M/(q+1));
        out = (n+1)*round(M/(q+1));
        theta = wrapTo2Pi((2*n-1)*atan2(-yy,-xx));
    end
    spiral((xx.^2+yy.^2)<=out.^2 & (xx.^2+yy.^2)>=in.^2) = exp(-1j*theta((xx.^2+yy.^2)<=out.^2 & (xx.^2+yy.^2)>=in.^2));
end

figure
imagesc(real(spiral));
colorbar

figure;
imagesc(abs(F(spiral)))

%% Analyzing Spiral Phase Mask as function of q
hi = Spiral_Phase_Mask(1,400);
figure
imagesc(real(hi))
