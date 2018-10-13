function Phase3D_G=MyMakingPhase_G(Nx,Ny,z,lambda,deltaX,deltaY,deltaZ)


k=1/lambda;

X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));

kx=repmat(X,1,Ny); % repeat X, 1 row Ny col
ky=repmat(Y,Nx,1); % repeat Y, Nx row 1 col
%kp=sqrt(kx.^2+ky.^2);
kp=kx.^2+ky.^2;

%term=k.^2-kp.^2;s
term=k.^2-kp;
term(term<0)=0;


%% Lei Tian commments:
% Phase3D_G = j*pi/lambda^2*exp(i*2*pi*z*kz)./kz 
% the extra pi constant is added for a redefinition of scattering
% potentional to be (n^2-n0^2)

 Phase3D_G = (1j*pi*(k^2)*deltaX)*exp(1j*2*pi*z*sqrt(term))./(sqrt(term));

end
