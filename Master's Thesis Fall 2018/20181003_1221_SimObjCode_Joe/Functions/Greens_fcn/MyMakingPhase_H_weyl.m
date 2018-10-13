function Phase3D_H = MyMakingPhase_H_weyl(Nx,Ny,z,lambda,deltaX,deltaY,...
    deltaZ,is_discr_fcn)

k=1/lambda;

X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));

kx=repmat(X,1,Ny); % repeat X, 1 row Ny col
ky=repmat(Y,Nx,1); % repeat Y, Nx row 1 col
kp=sqrt(kx.^2+ky.^2);
%kp=kx.^2+ky.^2;

term=k.^2-kp.^2;
%term=k.^2-kp;
term(term<0)=0;

% deltaZ is assumed to be equal to deltaX
Phase3D_H=(1i*pi*(k^2))*(exp(1i*2*pi*z*sqrt(term))./(sqrt(term)))*deltaX;

if (is_discr_fcn)
    % pixel transfer function
    P = sinc(deltaX*kx).*sinc(deltaY*ky);
    % the discretization transfer function
    Q = sinc(deltaX*kx).*sinc(deltaY*ky);%.*...
        %sinc(deltaz*((uu.^2+vv.^2)*lambda/2));
   
    Phase3D_H = Phase3D_H .* P .* Q; 
end
end