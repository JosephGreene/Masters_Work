function Phase3D_H = MyMakingPhase3D_H_weyl(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,~,is_rm_alias,is_discr_fcn)

if is_rm_alias
    Nx = Nx * 2;
    Ny = Ny * 2;
end

% axial plane z-coordinates
Z = [0:Nz-1].*deltaZ+offsetZ;
% camera is at the z=0 plane
Zcam = 0;

% propagation from all z-planes to camera plane at z=0 (i.e Zcam)
Phase3D_H=zeros(Nx,Ny,Nz);
for i=1:Nz
    Phase3D_H(:,:,i)=MyMakingPhase_H_weyl(Nx,Ny,abs(Zcam-Z(i)),lambda,deltaX,...
        deltaY,deltaZ,is_discr_fcn);
end