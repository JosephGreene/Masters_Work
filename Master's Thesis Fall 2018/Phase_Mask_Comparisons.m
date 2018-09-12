%Parameters for future use
%{
    GRIN lens is 1.8um across, smallest resolvable 'pixel' on mask will be
    ~4.7um. 
%}
%Progegation Equations
F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
prop = @(f,h) iF(F(f).*h);
lambda = 0.54; %Wavelength in um
dx = 4.5; %Real space 'pixel' size in um
M = 1800/4.5; %Number of pixels in diameter of GRIN lens
[xx,yy] = meshgrid([-M/2:M/2-1]*dx); %Spatial grid for use generating phase mask
[uu,vv] = meshgrid([-M/2:M/2-1]*(1/(M*dx))); %In Fourier plane, du is defined as 1/FOV

%Optimized coefficients utilized in generating phase masks in their
%respective papers.
%[1]
a = 37.53*10^-6;
b = 1.27;
%[2]
a2 = 6.5*10^-8;
%Proposed Phase Masks for extended depth of field uses
tangent = wrapTo2Pi(a*uu.^2*tan(b*xx)+a*yy.^2*tan(b*yy));
tangent((xx.^2+yy.^2) >= (0.95*M/2*dx).^2) = 0; %Mimic Aperture

cubic = a2*(xx.^3+yy.^3);
cubic((xx.^2+yy.^2) >= (0.95*M/2*dx).^2) = 0; %Mimic Aperture
wrapped = wrapTo2Pi(cubic);
binary = wrapped;
binary(binary >= pi) = 2*pi;
binary(binary < pi) = 0;

%Visualize masks
figure
imagesc(tangent)
title('Tangent Mask')
figure
imagesc(binary)
title('Binary Cubic Mask')

%% Generate PDF,OTF,MTF for on-axis cases
%{
    For a given pupil function p(u)
    iPSF @ z = hi(x,z) = |F{p(u).*exp(i*pi*lambda*z*u^2)}|^2 = h(x,z)^2
    OTF = F{h(x,z)} * F{h(x,z)} = autocorrelation of pupil function
    MTF = |OTF|
%}

%Tangent Phase Mask
OTF_tan = xcorr2(tangent,tangent);
MTF_tan = abs(OTF_tan);
figure
h = surf(OTF_tan);
set(h,'LineStyle','none');
figure
h = surf(MTF_tan);
set(h,'LineStyle','none');

%Binary Phase Mask
OTF_bin = xcorr2(binary,binary);
MTF_bin = abs(OTF_bin);
figure
h = surf(OTF_bin);
set(h,'LineStyle','none');
figure
h = surf(MTF_bin);
set(h,'LineStyle','none');

%% Defocus MTF of Tangential Phase Mask
defocus = [-6:2:6];
for i = 1:length(defocus)
    defocus_prop = exp(1i*pi*lambda*defocus(i)*(uu^2+vv^2));
    %Tangent Phase Mask
    OTF_tan = xcorr2(tangent.*defocus_prop,tangent.*defocus_prop);
    MTF_tan = abs(OTF_tan);
    figure
    h = surf(MTF_tan);
    set(h,'LineStyle','none');
    title(strcat('Tan Mask, defocus =',int2str(defocus(i))))
    
    %Binary Phase Mask
    OTF_bin = xcorr2(binary.*defocus_prop,binary.*defocus_prop);
    MTF_bin = abs(OTF_bin);
    figure
    h = surf(MTF_bin);
    set(h,'LineStyle','none');
    title(strcat('Binary Mask, defocus =',int2str(defocus(i))))
    
end
%% References
%{
[1] Optimized asymmetrical tangent phase mask to obtain defocus invariant modulation transfer function in incoherent imaging systems
[2] A novel two- and multi-level binary phase mask designed for enhanced depth-of-focus
%}