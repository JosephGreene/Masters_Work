function grad=MyGradient(resid,f,E_rb_all,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_order,F_nopad,F_pad,Ft,is_rm_alias,sv_sq_H,sv_sq_G)

if is_rm_alias
    F = F_pad;
else
    F = F_nopad;
end

f_temp = f;
% Compute gK = [Hh r] .* conj(uK)
%% Lei Tian comments:
% 1) the physical meaning of this step is to "back-propagate" the residual to 
%    the 3D object space, thus conj(.) is used.
% 2) the adjoint needs take into account the 2x factor in the 
%    hologram formation equation 2*real(E).
%    essentially, this is a quasi-Newton's method by taking into account
%    the non-unit svd values in the Forward model and try to fix it by
%    parts
cEsp=F(resid);

cEs=conj(phase3D_H).*repmat(cEsp,[1 1 Nz]) * 1/1;

% eta here denotes the virtual source due to the residual
eta=zeros(Nx,Ny,Nz);
for i=1:Nz
    tmp = ifft2(ifftshift(cEs(:,:,i)));
    eta(:,:,i)= tmp(1:Nx,1:Ny);
end

gK = conj(E_rb_all(:,:,:,born_order)).*eta; 
% scaling with singular value sq
gK = gK/sv_sq_H;
% if born_order=3, we have K=2,g2, & E_rb corresponding to u2


if(born_order > 1)
    
    % Compute vK = [Hh r] .* conj(f)
    vK = conj(f_temp).*eta;
    % scaling with singular value sq
    vK = vK/sv_sq_H;
    
    g_k_plus_1 = gK; %gK=g2
    v_k_plus_1 = vK;
    
    k = born_order;
    while (k > 1) %k=3 %k=2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % % Compute [Gh vk+1] .* conj(uk)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        cEsp_3D = zeros(2*Nx,2*Ny,Nz);
        for i = 1:Nz
            cEsp_3D(:,:,i)=fftshift(fft2(v_k_plus_1(:,:,i),2*...
                size(v_k_plus_1(:,:,i),1),2*size(v_k_plus_1(:,:,i),2)));
            % cEsp_3D(:,:,i) = F(v_k_plus_1(:,:,i));
        end

        cEs_4D = zeros(2*Nx,2*Ny,Nz,Nz);
        for i = 1:Nz
            cEs_4D(:,:,:,i)=conj(phase3D_G(:,:,:,i))...
                .*repmat(cEsp_3D(:,:,i),[1 1 Nz]) * 1/1;
        end
        
        cEs_sum = sum(cEs_4D,4);

        eta=zeros(Nx,Ny,Nz);
        for i=1:Nz
            tmp = ifft2(ifftshift(cEs_sum(:,:,i)));
            eta(:,:,i) = tmp(1:Nx,1:Ny);
            %eta(:,:,i)=Ft(cEs_sum(:,:,i));
        end

        %% Lei Tian commments:
        % need to consider the case when abs(E_rb_all) is not 1?
        %    Note: essentially, this is a quasi-Newton's method by taking into account
        %    the non-unit svd values in the Forward model and try to fix it by
        %    parts
        G_v_k_plus_1=eta.*conj(E_rb_all(:,:,:,k-1)); 
        % scaling with singular value sq
        G_v_k_plus_1 = G_v_k_plus_1/sv_sq_G;
        
        %born_order=3, we have k=3 for 1st iter, k-1 = 2, corresponding to u1  
        % E_rb corresponds to k-1=1 i.e. u0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Compute gk = gk+1 + [Gh vk+1] .* conj(uk)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        gk = g_k_plus_1 + G_v_k_plus_1; %g1 %g0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Compute vk = [Gh vk+1] .* f
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vk = conj(f_temp).*eta; %v1 v0
        % scaling with singular value sq
        vk = vk/sv_sq_G;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Recursion of g & v variablezs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        g_k_plus_1 = gk;
        v_k_plus_1 = vk;    

        k = k - 1; %k = 2
    end
else
    gk = gK;
end

grad = gk;

end