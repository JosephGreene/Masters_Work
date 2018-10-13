function E=MyFieldsBackPropagation(E0,Nx,Ny,Nz,Phase3D_Prop,F,Ft)

E=zeros(Nx,Ny,Nz);
% Compute fft2 of E0 [which is ones(Nx,Ny)]
cE0=F(E0);

for i=1:Nz
    % Multuply fft2 of E0 by conjugte of [Phase3D (complex exponential) * Pupil (ones(Nx,Ny)]
    % Convolve with Phase3D in spatial domain
    cE=cE0.*Phase3D_Prop(:,:,i); % back prop, conj
    
    % Compute E as IFT of above multiplication in fourier domain
    %(convolution in spatial domain)
    E(:,:,i)=Ft(cE);
end