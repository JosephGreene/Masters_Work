%Parameters for future use
%{
    GRIN lens is 1.8um across, smallest resolvable 'pixel' on mask will be
    ~4.7um. 

    ****BU miniscope pixel size in relity is 6um but this causes undersampling
    so I changed it in this system****
%}
%% Basic Parameters
%Choose Which Masks to generate
showTan = 0;
showCub = 0;
showAxi = 1;

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));
prop = @(f,h) iF(F(f).*h);

%Physical parameters of miniscope (unchagable) in real space
NA0 = 0.45; %NA of the system -> range of ATF (2*NA is range of OTF)
f = 3.3; %mm
dimag = 10; %mm
lambda = 0.00054; %Wavelength in mm
Nyquest = lambda/(4*NA0); %Nyquest requirement for incoherent system for comparison
pixel = 0.003; %ideal pixel size in mm
mag = 10; %magnification of system
dx = pixel/mag; %Real space 'pixel' size in mm
N = 400; %desired # of pixels
fov = N*dx; %N = 1000 in this case
[xx,yy] = meshgrid([-N/2:N/2-1]*dx); %Spatial grid for use generating phase mask
DoF = (lambda/(2*NA0^2)); %Native depth of field in this system
%% Fourier space
du = 1/fov;
[uu,vv] = meshgrid([-N/2:N/2-1]*du); %In Fourier plane, du is defined as 1/FOV
NAx = uu*lambda; %converting to NA space (alottable angles) -> unitless! easy of scaling and design
NAy = vv*lambda;
NA = sqrt(NAx.^2+NAy.^2); %NA @ any given point
%% Generate Masks as function of NA (unitless for optimal flexibility)

%Optimized coefficients utilized in generating phase masks in their
%respective papers.
%[1] tangent mask
a = 37.53*10^-2;
b = 1.27;
%[2] Binarized Cubic Phase Mask
a2 = 6.5*10^2;
%[3] Binarized Axicon
% alpha < NA^2*f/lambda converges spherical waves after lens to generate
% bessel-like beam
alpha = NA0^2*f/(3*lambda)*10^-2;

%Proposed Phase Masks for extended depth of field uses
if showTan
    tangent = wrapTo2Pi(a*NAx.^2*tan(b*NAx)+a*yy.^2*tan(b*NAy));
    tangent((NAx.^2+NAy.^2) >= (0.45)^2) = 0; %Mimic Aperture
    %Visualize masks
    figure
    imagesc(tangent)
    title('Tangent Mask')
end

if showCub
    cubic = a2*(NAx.^3+NAy.^3);
    cubic((NAx.^2+NAy.^2) >= (0.45)^2) = 0; %Mimic Aperture
    wrapped = wrapTo2Pi(cubic);
    bincubic = wrapped;
    bincubic(bincubic >= pi) = 2*pi;
    bincubic(bincubic < pi) = 0;
    figure
    imagesc(bincubic)
    title('Binary Cubic Mask')
end

if showAxi
    axicon = exp(-2*pi*1i*alpha*sqrt((NAx.^2+NAy.^2))); %Paper suggests normalizing pupil coordinate
    binaxi = axicon;
    binaxi(real(binaxi) < 0) = -1;
    binaxi(real(binaxi) >= 0) = 1;
    binaxi((NAx.^2+NAy.^2) >= (0.45)^2) = 0;
    figure
    imagesc(binaxi)
    title('Binary Axicon Mask')
end

%% Generate PDF,OTF,MTF for on-axis cases
%{
    For a given pupil function p(u)
    iPSF @ z = hi(x,z) = |F{p(u).*exp(i*pi*lambda*z*u^2)}|^2 = h(x,z)^2
    OTF = F{h(x,z)} * F{h(x,z)} = autocorrelation of pupil function
    MTF = |OTF|
%}
if showTan
    %Tangent Phase Mask
    OTF_tan = xcorr2(tangent,tangent);
    MTF_tan = abs(OTF_tan);
    figure
    h = surf(OTF_tan);
    set(h,'LineStyle','none');
    figure
    h = surf(MTF_tan);
    set(h,'LineStyle','none');
end

if showCub
    %Binary Phase Mask
    OTF_cub = xcorr2(bincubic,bincubic);
    MTF_cub = abs(OTF_cub);
    figure
    h = surf(OTF_cub);
    set(h,'LineStyle','none');
    figure
    h = surf(MTF_cub);
    set(h,'LineStyle','none');
end

if showAxi
    %Binary Axicon
    OTF_axi = xcorr2(binaxi,binaxi);
    MTF_axi = abs(OTF_axi);
    figure
    h = surf(OTF_axi);
    set(h,'LineStyle','none');
    figure
    h = surf(MTF_axi);
    set(h,'LineStyle','none');
end
%% Defocus MTF of Tangential Phase Mask
%{
    For the miniscope, the focal length in 3.3mm, the distance to the
    imager is 10mm so the in focus distance for the object is 4.92 or ~
    5mm. So our dz defocus should be centered @ 5
%}
dz = [-0.02:0.0025:0.02];
MTF_in_focus = ones(2*N-1);
ip_lens = zeros(1,length(dz));
ip_tan = zeros(1,length(dz)); %Preallocate variable to hold inner product
ip_cub = zeros(1,length(dz));
ip_axi = zeros(1,length(dz));

for i = 1:length(dz)
    defocus = dz(i);
    defocus_prop = exp(1i*pi*lambda*defocus.*(uu.^2+vv.^2));
    
    %imagesc(real(defocus_prop))
    %break
    if showTan
        %Tangent Phase Mask: Calculate OFT/MTF with defocus
        %dOTF_tan = F(iF(tangent.*defocus_prop).*iF(tangent.*defocus_prop));
        dOTF_tan = xcorr2(tangent.*defocus_prop,tangent.*defocus_prop);
        dMTF_tan = abs(dOTF_tan);
%         figure
%         h = surf(MTF_tan);
%         set(h,'LineStyle','none');
%         title(strcat('Tan Mask, defocus =',num2str(dz(i))))
        %Calculate frobenius inner product of in focus and out of focus MTF
        ip_tan(i) = trace(MTF_tan.*transpose(dMTF_tan));
    end
    
    if showCub
        %Binary cubic Phase Mask
        dOTF_cub = xcorr2(bincubic.*defocus_prop,bincubic.*defocus_prop);
        dMTF_cub = abs(dOTF_cub);
%         figure
%         h = surf(MTF_bin);
%         set(h,'LineStyle','none');
%         title(strcat('Binary Mask, defocus =',num2str(dz(i))))
        ip_cub(i) = trace(MTF_cub.*transpose(dMTF_cub));
    end
    
       if showAxi
        %Binary Axicon Phase Mask
        dOTF_axi = xcorr2(binaxi.*defocus_prop,binaxi.*defocus_prop);
        dMTF_axi = abs(dOTF_axi);
%         figure
%         h = surf(MTF_bin);
%         set(h,'LineStyle','none');
%         title(strcat('Binary Mask, defocus =',num2str(dz(i))))
        ip_axi(i) = trace(MTF_axi.*transpose(dMTF_axi));
       end
       
       
       
       %plot(1:length(dMTF_axi),dMTF_axi(0,:))
       %legend(strcat('Axi, defocus =',num2str(dz(i))))
       
end

figure
hold on
plot(ip_tan./max(ip_tan),'r')
plot(ip_cub./max(ip_cub),'b')
plot(ip_axi./max(ip_axi),'k')
hold off
legend('tan','binary cubic','binary axicon')
%% References
%{
[1] Optimized asymmetrical tangent phase mask to obtain defocus invariant modulation transfer function in incoherent imaging systems
[2] A novel two- and multi-level binary phase mask designed for enhanced depth-of-focus
[3] Design of binary phase filters for depth-of-focus extension via
binarization of axisymmetric aberrations
%}