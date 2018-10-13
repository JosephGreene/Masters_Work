function Phase3D_Prop=MyMakingPhase_Prop(Nx,Ny,z,lambda,deltaX,deltaY)

k=1/lambda;

X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));

kx=repmat(X,1,Ny); % repeat X, 1 row Ny col
ky=repmat(Y,Nx,1); % repeat Y, Nx row 1 col
%kp=sqrt(kx.^2+ky.^2);
kp=kx.^2+ky.^2;

%term=k.^2-kp.^2;
term=k.^2-kp;
term(term<0)=0;

%% Lei Tian Comments:
% Angular spectrum representation of forward propagator, which is different
% from the Green's function angular spectrum expansion (also known as 
Phase3D_Prop=exp(1j*2*pi*z*sqrt(term));
end