classdef FwdInvModel < handle
    %%% A - Forward operator of the form:
    %%% z = A(x)
    %%% where, z is the estimated hologram
    %%%        x is the 3D object scattering density (n^2-1)
    %%%
    %%% AT_init - Transpose operator of the form
    %%% x = AT(y)
    %%% where, y is the hologram measurement
    %%%        x is the least sq solution of the object
    %%%
    %%% Grad - Gradient operator of the form
    %%% grad = Grad(resid)
    %%% where, resid is the residual(y-z)
    %%%        grad is the gradient of data fidelity D=|y-z(x)|^2 wrt x
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sigSize; % size of the signal
        measSize; % size of the measurements
        A; % forward operator
        Grad; % inverse operator
        AT_init; % initialization operator
        step; % step size for gradient algorithms
        bornorder;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = FwdInvModel(A,Grad,AT_init,bornorder,sigSize, sv_sq)
            %%% Constructor
            % Set parameters
            this.A = A;
            this.Grad = Grad;
            this.AT_init = AT_init;
            this.sigSize = sigSize;
            this.bornorder = bornorder;
            % Step-size is 1/max(eig(A'A)) or 1/max(singular_val(A))^2
            if(bornorder > 1)
                this.step = 1/sv_sq * 30;
            else
                this.step = 1/sv_sq * 30;
            end
        end
    end
end