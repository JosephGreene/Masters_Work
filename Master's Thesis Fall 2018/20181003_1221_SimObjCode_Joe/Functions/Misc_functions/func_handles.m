% Defining FFT function handles
F = @(x) fftshift(fft2(x));
Fpad = @(x) fftshift(fft2(x,2*size(x,1),2*size(x,2)));
Ft = @(x) ifft2(ifftshift(x));