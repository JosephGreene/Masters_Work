%xhatnext             = proxTV(xhatnext, step*tau , gradObj);
function [xhat, outs] = proxTV(y       , lam      , gradObj, varargin)

%%% TV denoising of a complex image 
%%%
%%% U. S. Kamilov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeSnr = @(x, xhat) 20*log10(norm(MyV2C(x))/norm(MyV2C(x)-MyV2C(xhat)));
tvCost = @(x) sum(sum(sum(sqrt(sum(abs(gradObj.mult(x)).^2, 4)))));
computeCost = @(x) 0.5*(norm(MyV2C(x)-MyV2C(y))^2)+lam*tvCost(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numIter = 10; % Number of iterations
plotRecon = false; % Plot the progress
tol = 1e-6; % Tolerance for relative change
verbose = false; % Display messages

rho = 1; % quadratic parameter
xhat0 = y; % initial image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = length(varargin); % Number of options

for i = 1:2:nargs % Go through options
    % Extract name/value pair
    name = lower(varargin{i});
    value = varargin{i+1};
    
    switch(name)
        
        case 'numiter'
            numIter = value;
        case 'plotrecon'
            plotRecon = value;
        case 'rho'
            rho = value;
        case 'tol'
            tol = value;
        case 'x'
            x = value;
        case 'xhat0'
            xhat0 = value;
        case 'verbose'
            verbose = value;
        otherwise
            error('denTvAdmm: input is not recognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isTrueXSet = exist('x', 'var'); % oracle

% initial values
xhat = xhat0;
d = gradObj.mult(xhat);
s = zeros(size(d));

% initialize arrays holding per-iteration cost and snr
outs.cost = zeros(numIter, 1);
if(isTrueXSet)
    outs.snr = zeros(numIter, 1);
end

% create the figure
if(plotRecon)
    figHandle = figure('Color', 'w', 'Name', 'denTvAdmm');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIter = 1:numIter
    
    % message to command line
    if(verbose)
        fprintf('[proxTV: %d/%d]', iIter, numIter);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update estimates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xhatprev = xhat; % store the current estimate before update
    
    % update d
    d = shrink3D(gradObj.mult(xhat)-s/rho, lam/rho);
    
    % update xhat
    dat = y + rho*gradObj.multTranspose(d + s/rho);
    xhat = ifftn(fftn(dat) ./ (1 + rho*gradObj.freqResponseDTD));
    
    % update s
    s = s + rho*(d - gradObj.mult(xhat));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Illustrate results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(isTrueXSet)
        outs.snr(iIter) = computeSnr(x, xhat);
    end
    outs.cost(iIter) = computeCost(xhat);
    
    % print to command line
    if(verbose)
        fprintf('[cost: %.1e]', outs.cost(iIter));
        if(isTrueXSet)
            fprintf('[SNR: %.2f]', outs.snr(iIter));
        end
        fprintf('\n');
    end
    
    % plot
    if(plotRecon)
        figure(figHandle);
        set(figHandle,...
            'Name', sprintf('[proxTV: %d/%d]', iIter, numIter));
        
        if(isTrueXSet)
            subplot(1, 2, 1);
            semilogy(1:iIter, outs.cost(1:iIter), '-b',...
                iIter, outs.cost(iIter), 'ro', 'LineWidth', 1.2);
            title(sprintf('cost: %.1e', outs.cost(iIter)));
            xlim([1 numIter]);
            grid on;
            set(gca, 'FontSize', 14);
            
            subplot(1, 2, 2);
            plot(1:iIter, outs.snr(1:iIter), '-b',...
                iIter, outs.snr(iIter), 'ro', 'LineWidth', 1.2);
            title(sprintf('snr: %.2f dB', outs.snr(iIter)));
            xlim([1 numIter]);
            grid on;
            set(gca, 'FontSize', 14);
        else
            semilogy(1:iIter, outs.cost(1:iIter), '-b',...
                iIter, outs.cost(iIter), 'ro', 'LineWidth', 1.2);
            title(sprintf('cost: %.1e', outs.cost(iIter)));
            xlim([1 numIter]);
            grid on;
            set(gca, 'FontSize', 14);
        end
        drawnow;
    end
    
    % check tolerance
    if(norm(xhat(:)-xhatprev(:))/norm(xhatprev(:)) < tol)
        outs.cost = outs.cost(1:iIter);
        if(isTrueXSet)
            outs.snr = outs.snr(1:iIter);
        end
        break;
    end
end