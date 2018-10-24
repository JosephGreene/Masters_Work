function [n1,n2] = Signc2(N)
%{
Produces 2D meshgrid of discretized, signum corrdinates as described in
"An approach to the 2D Hilbert Transform for Image Processing 
Applications". Assumes image is NxN. Signum-Combined-2D
%}
[xx,yy] = meshgrid(0:N-1);
xx(xx <= floor(N/2)-1) = 1;
xx(xx > floor(N/2)-1 & xx <= N-1) = -1;
xx(:,1) = 0;
xx(:,N) = 0;

yy(yy <= floor(N/2)-1) = 1;
yy(yy > floor(N/2)-1 & xx <= N-1) = -1;
yy(1,:) = 0;
yy(N,:) = 0;

n1 = xx;
n2 = yy;
end