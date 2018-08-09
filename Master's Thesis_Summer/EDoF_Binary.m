%% Extended Depth of Field Binary Phase Masks
%% Generalized Cubic Mask

N = 400;
[xx,yy] = meshgrid(-N/2:N/2);
cubic = 6.5*(xx.^3+yy.^3)*10^-6;
cubic((xx.^2+yy.^2) >= (0.95*N/2).^2) = 0; %Mimic Aperture
wrapped = wrapTo2Pi(cubic);

binary = wrapped;
binary(binary >= pi) = 2*pi;
binary(binary < pi) = 0;

figure
imagesc(cubic)
title('Cubic Phase Mask')
colorbar

figure
imagesc(wrapped)
title('Phase Wrapped Cubic')
colorbar

figure
imagesc(binary)
title('Binary Wrapped Cubic')
colorbar
