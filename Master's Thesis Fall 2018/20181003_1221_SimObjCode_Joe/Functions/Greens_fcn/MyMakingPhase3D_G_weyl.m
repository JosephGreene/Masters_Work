function Phase3D_G = MyMakingPhase3D_G_weyl(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,~,is_rm_alias,~)

if is_rm_alias
    Nx = Nx * 2;
    Ny = Ny * 2;
end

% axial plane z-coordinates
Z = [0:Nz-1].*deltaZ+offsetZ;

Phase3D_G = zeros(Nx,Ny,Nz,Nz);
for j=1:Nz % separate Phase3D for each z
    for i=1:Nz
        if i~=j
            % propagation from slice at z(i) to sclice at z(j)
            Phase3D_G(:,:,i,j)=...
                MyMakingPhase_G(Nx,Ny,abs(Z(j)-Z(i)),lambda,deltaX,deltaY,deltaZ);
        else
            Phase3D_G(:,:,i,j)=zeros(Nx,Ny);
        end
            
    end
end