function  [Phase3D_H, Greens_spatial, dist] = MyMakingPhase3D_H_conv(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,is_discr_fcn)

if is_rm_alias
    pad_len = 2;
    Nxs = Nx * 2;
    Nys = Ny * 2;
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
end

k=(2*pi)/lambda;

% x-y-z axes coordinates
X = [ceil(-Nx/2):1:ceil(Nx/2-1)]'*deltaX;
Y = [ceil(-Ny/2):1:ceil(Ny/2-1)]*deltaY;
Z = [0:Nz-1].*deltaZ+offsetZ;

X_mat=repmat(X,1,Ny); % repeat X, 1 row Ny col
Y_mat=repmat(Y,Nx,1); % repeat Y, Nx row 1 col

% calculating phase shift in x and y
kx = 0:Nx-1;
ky = 0:Ny-1;
phase_shiftx = exp(+1i*(2*pi*kx/(Nx))*(Nx/2)); 
phase_shifty = exp(-1i*(2*pi*ky/(Ny))*(Ny/2));
Phase_x = repmat(phase_shiftx,Nx,1); % FFT compensated for phase shift in  x
Phase_y = repmat(phase_shifty',1,Ny); % FFT compensated for phase shift in  y

% propagation from all z-planes to camera plane at z=0 (i.e Zcam)

Phase3D_H=zeros(Nx,Ny,Nz);
Greens_spatial = zeros(Nx,Ny,Nz);
dist = zeros(Nx,Ny,Nz);
debg = 1; % debugging 
for i=1:Nz
    Z_mat = repmat(Z(i),Nx,Ny);
    term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
    if(debg)
        dist(:,:,i) = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
        Greens_spatial(:,:,i) = exp(1i*k*term)./term;
    end
    % shifting because matlab starts counting from 0, whereas, we want
    % the counting to start form -ve. Multiplying with const term k^2/4pi
    % for object scattering density and dV for the area under unit impulse
    const = ((k^2)/(4*pi))*deltaX^3;
    Phase3D_H(:,:,i) = (F(exp(1i*k*term)./term).*Phase_x.*Phase_y)*const;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (is_discr_fcn)
    Xs=[ceil(-Nxs/2):1:ceil(Nxs/2-1)]'.*(1/(Nxs*deltaX));
    Ys=[ceil(-Nys/2):1:ceil(Nys/2-1)].*(1/(Nys*deltaY));
    kxs=repmat(Xs,1,Nys); % repeat X, 1 row Ny col
    kys=repmat(Ys,Nxs,1); % repeat Y, Nx row 1 col
    % pixel transfer function
    P = sinc(deltaX*kxs).*sinc(deltaY*kys);
    % the discretization transfer function
    Q = sinc(deltaX*kxs).*sinc(deltaY*kys);%.*...
        %sinc(deltaZ*((kxs.^2+kys.^2)*lambda/2));
   
    Phase3D_H = Phase3D_H .* repmat(P,[1 1 Nz]) .* repmat(Q,[1 1 Nz]); 
end

