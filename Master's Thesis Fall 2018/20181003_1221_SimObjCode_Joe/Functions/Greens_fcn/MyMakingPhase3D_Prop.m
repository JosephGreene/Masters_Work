function Phase3D_Prop = MyMakingPhase3D_Prop(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ)

% axial plane z-coordinates
Z = [0:Nz-1].*deltaZ+offsetZ;
% camera is at the z=0 plane
Zcam = 0;

Phase3D_Prop=zeros(Nx,Ny,Nz);
for i=1:Nz
    Phase3D_Prop(:,:,i)=...
        MyMakingPhase_Prop(Nx,Ny,Zcam-Z(i),lambda,deltaX,deltaY);
end
end