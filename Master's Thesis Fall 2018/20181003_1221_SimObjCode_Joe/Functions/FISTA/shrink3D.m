function xhat = shrink3D(y, tau)
%%% Shrinkage function for complex 3D-TV denoising.
%%% U. S. Kamilov

% shortcut for computing norm along the fourth dimension
computeNormInZ = @(x) sqrt(sum(abs(x).^2, 4));

% convert data to a 6D array with real elements
yv = zeros([size(y, 1), size(y, 2), size(y, 3), 6]);
yv(:,:,:,1) = real(y(:,:,:,1));
yv(:,:,:,2) = imag(y(:,:,:,1));
yv(:,:,:,3) = real(y(:,:,:,2));
yv(:,:,:,4) = imag(y(:,:,:,2));
yv(:,:,:,5) = real(y(:,:,:,3));
yv(:,:,:,6) = imag(y(:,:,:,3));

% compute the norm in dimension 3
norm_yv = computeNormInZ(yv);

amp = max(norm_yv-tau, 0);
norm_yv(norm_yv<=0) = 1; % to avoid division by 0

xhatv = repmat(amp./norm_yv, [1, 1, 1, 6]) .* yv;

xhat = zeros(size(y));
xhat(:,:,:,1) = xhatv(:,:,:,1) + 1j*xhatv(:,:,:,2);
xhat(:,:,:,2) = xhatv(:,:,:,3) + 1j*xhatv(:,:,:,4);
xhat(:,:,:,3) = xhatv(:,:,:,5) + 1j*xhatv(:,:,:,6);