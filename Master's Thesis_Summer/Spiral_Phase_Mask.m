function spiral = Spiral_Phase_Mask(n,N)
%{
Spiral Phase Mask described in "Spiral Phase Contrast Imaging in Microscopy"
%}
spiral = zeros(N,N);
M = round(N/2);
[xx,yy] = meshgrid([-M:M-1]);
theta = atan2(yy,xx);

for q = 1:n
    if(q == 1)
        mask = (xx.^2+yy.^2) <= (2/(n+1)*(N-M-1))^2;
        spiral(mask) = exp(1j*theta(mask));
    else
        mask = (xx.^2+yy.^2) <= ((q+1)/(n+1)*(N-M+1))^2 & (xx.^2+yy.^2) >= ((q)/(n+1)*(N-M+1))^2;
        spiral(mask) = exp(1j*wrapTo2Pi(q*theta(mask)));
    end
    spiral(M,M) = 0; %As described by paper
end