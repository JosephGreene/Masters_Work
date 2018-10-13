clear all; close all; clc; add_paths; func_handles;
% Compressive Holography with Recursive Born
% Solver: FISTA (Ulugbek's FPPA)


% Set Parameters215
super_res_factor = 1; % super resolution in object domain w.r.t pixel size
is_rm_alias = 1;      % remove aliasing in convolutions? currenlty the code
%only works when this is set to 1. Hardcoding is in the parfor loops
is_discr_fcn = 0;     % LPF due to object domain discretization
is_perf_metric = 0;

doBPM = 1;          % flag for backpropgation (least sq) result
max_born_order = 1; % max born order for computation

%%% Main Computation Loop
for tau = [1]
    for born_order = [1]

        tic

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Hologram preprocessing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load data (matfile)
        load('img(1)_filt.mat');
        I = I_nobg;
        clear I_nobg;
        % Load image file (only for display)
        %I_u16 = imread('img(1)_filt.tif');
        
        % crop 256x256 region (computation for starting pt on paper)
        holo_size = 256;
        pt = (2048-holo_size)/2;
        Holo_meas = I(pt:pt+(holo_size-1),pt:pt+(holo_size-1));
        %Holo_meas_disp = I_u16(pt:pt+(holo_size-1),pt:pt+(holo_size-1));
        clear I I_u16;
        
        % no upsampling using upsample4x(Holo_meas)
        
        % DC subtraction
        Avg_int = mean2(Holo_meas);
        Avg_amp = sqrt(Avg_int);
        Holo = (Holo_meas-Avg_int)/(2*Avg_amp);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initializing simulation parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nx, ny]=size(Holo_meas);   % no. of x-y pixels
        sim_params_5_5;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialize incident field and propagation operators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % define incident field magnitude
        % Lei Tian edits:
        % mean2(Holo) is the estimated DC intensity of input field
        E0  = ones(nx,ny).*Avg_amp;
        
        % Phase3D_H defines propagation from 3D object to 2D hologram
        Phase3D_H = MyMakingPhase3D_H_weyl(nx,ny,nz,lambda,...
            deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,is_discr_fcn);
        
        % Phase3D_G defines propagation from 3D object to within the object
        % (repetitive convolution with this yeilds higher order fields)
        Phase3D_G = MyMakingPhase3D_G_weyl(nx,ny,nz,lambda,...
            deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,is_discr_fcn);
        
        % Phase3D_Prop defines unhindered propagation of E0 within the object
        % Lei Tian comments:
        % field propagatino needs to be done based on Angular spectrum,
        % which is slightly different from the Green's function.
        Phase3D_Prop = MyMakingPhase3D_Prop(nx,ny,nz,lambda,...
            deltaX,deltaY,deltaZ,offsetZ);

        % back-propagating E0 in free space from z=0 to within obj cube
        E = MyFieldsBackPropagation(E0,nx,ny,nz,Phase3D_Prop,F,Ft);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function handles for forward (A) and inverse (AT) oprators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % forward operator
        A = @(f) MyForwardPropagation_rb(f,E,nx,ny,nz,Phase3D_H,...
        Phase3D_G,born_order,F,Fpad,Ft,is_rm_alias,is_perf_metric,deltaZ,super_res_factor,deltaX);

        % compute singular values of A
        % singular values of A are sqrt(eig(A'A))
        % where, eig(A'A) = 2*(Avg_amp*nz*pi^2*deltaZ^2)/lambda^2
        sv_sq = (Avg_int*nz*pi^2*deltaZ^2)/lambda^2;

        
        % adjoint/backward operator
        Grad = @(resid,f,E_rb_all) MyGradient(resid,f,E_rb_all,nx,ny,nz,...
            Phase3D_H,Phase3D_G,born_order,F,Fpad,Ft,is_rm_alias,1,1);

        % AT_init does not depend on f, AT does. Used for initial guess of f
        AT_init = @(g) MyAdjointPropagation(g,E,nx,ny,nz,Phase3D_H,Fpad,...
            Ft,sv_sq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % backpropagation reconstruction (least sq solution)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        results_dir = sprintf('Results/tau_%1.0e',tau);
        mkdir(results_dir);
        fileName = sprintf('%s/prox_descent_data_%1.0e_b%d.txt',results_dir,tau,born_order);
        fileID = fopen(fileName,'wt');
        
        doBPM = 0;
        if doBPM
            % computes backpropagation (only in the 1st iteration)
            fprintf('Backpropagation...\n');
            fprintf(fileID,'Backpropagation...\n');
            
            % compute least sq solution by backpropagation
            BPField = AT_init(Holo);
            %BPField = reshape(MyV2C(f_least_sq),nx,ny,nz);
            
            % Create folders for saving BPM results
            bpm_dir = 'Results/BPM';
            mkdir(bpm_dir);
            mkdir([bpm_dir,'/BPM_abs']);
            mkdir([bpm_dir,'/BPM_real']);
            mkdir([bpm_dir,'/BPM_imag']);
            
            % BPM Abs
            fn =[bpm_dir,'/BPM_abs/','BPM_abs2_'];
            bp_min= min(abs(BPField(:)));
            bp_max = max(abs(BPField(:)));
            for m =1 : nz
                slice = mat2gray(abs(BPField(:,:,m)),[bp_min bp_max]);
                filename = [fn, num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end
            
            % BPM Real
            fn =[bpm_dir,'/BPM_real/','BPM_real2_'];
            bp_min= min(real(BPField(:)));
            bp_max = max(real(BPField(:)));
            for m =1 : nz
                slice = mat2gray(real(BPField(:,:,m)),[bp_min bp_max]);
                filename = [fn, num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end
            
            % BPM Imag
            fn =[bpm_dir,'/BPM_imag/','BPM_imag2_'];
            bp_min= min(imag(BPField(:)));
            bp_max = max(imag(BPField(:)));
            
            % Save BPM reconstruction as bmp images
            for m =1 : nz
                slice = mat2gray(imag(BPField(:,:,m)),[bp_min bp_max]);
                filename = [fn, num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end
            % no BPM for next iters
            doBPM = 0; 
            % clean up memory
            clear BPField E0 Phase3D_Prop
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % Run FISTA to get final result. Ulugbek's implementation 'FPPA'.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Generate forward model
        fwd_inv_obj = FwdInvModel(A,Grad,AT_init,born_order,[nx,ny,nz],sv_sq);

        %%% Parameter selection
        numIter = 150; 
        
        %%% Run FISTA
        %fprintf('compressive reconstruction with FPPA\n\n');
        fprintf(fileID,'compressive reconstruction with FPPA\n\n');
        [f, ~] = fistaEst(fwd_inv_obj, tau, Holo,fileID, results_dir);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Saving FISTA Results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % make output directories for results
        f_outdir = [results_dir,'/born',num2str(born_order),'/iters_all'];
        mkdir(f_outdir);
        mkdir([f_outdir,'/f_abs']);
        mkdir([f_outdir,'/f_real']);
        mkdir([f_outdir,'/f_imag']);

        %%% f-reconstruction result - saving
        if 1 
            % Save f as mat file
            fn =[f_outdir,'/f_born',num2str(born_order),'.mat'];
            save(fn,'f')

            % Save abs(f) as image files
            fn =[f_outdir,'/f_abs'];
            f_min= min(abs(f(:)));
            f_max = max(abs(f(:)));
            for m =1 : nz
                slice = mat2gray(abs(f(:,:,m)),[f_min,f_max]);
                filename = [fn,'/TV_born',num2str(born_order),'_abs_', num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end

            % Save real(f) as image files
            fn =[f_outdir,'/f_real'];
            f_min= min(real(f(:)));
            f_max = max(real(f(:)));
            for m =1 : nz
                slice = mat2gray(real(f(:,:,m)),[f_min,f_max]);
                filename = [fn,'/TV_born',num2str(born_order),'_real_', num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end

            % Save imag(f) as image files
            fn =[f_outdir,'/f_imag'];
            f_min= min(imag(f(:)));
            f_max = max(imag(f(:)));
            for m =1 : nz
                slice = mat2gray(imag(f(:,:,m)),[f_min,f_max]);
                filename = [fn,'/TV_born',num2str(born_order),'_imag_', num2str(m),'.bmp'];
                imwrite(slice,filename,'bmp')
            end
        end

        %%% Measured hologram - saving
        figholo = 0;
        if figholo
            figure;
            % Waleed: show captured 2D hologram
            gca = imagesc(abs(Holo_meas));
            title('Measured Hologram');
            colorbar;
            filename =[f_outdir,'/','Measured Hologram.png'];
            saveas(gcf,filename,'png')
        end
        close all;

        % print time in minutes
        toc/60
    end
end

