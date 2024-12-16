function [xout, bestval] = test_DS(example_idx,nsamples,d,Lz)
%test for differential evolution
objFctHandle = @(Params,xs)fxval(xs,d,example_idx);

% Define parameter names, ranges and quantization:

% 1. column: parameter names
% 2. column: parameter ranges
% 3. column: parameter quantizations
% 4. column: initial values (optional)

paramDefCell = {
	'parameter1', [-Lz Lz], 0.01
	'parameter2', [-Lz Lz], 0.01
    'parameter3', [-Lz Lz], 0.01
    'parameter4', [-Lz Lz], 0.01
};

% Set initial parameter values in struct objFctParams 
objFctParams.parameter1 = 4; objFctParams.parameter2 = 4;
objFctParams.parameter3 = 4; objFctParams.parameter4 = 4;

% Set single additional function parameter
objFctSettings = 100;

% Get default DE parameters
DEParams = getdefaultparams;

% Set number of population members (often 10*D is suggested) 
DEParams.NP = 10*d;

% Do not use slave processes here. If you want to, set feedSlaveProc to 1 and
% run startmulticoreslave.m in at least one additional Matlab session.
DEParams.feedSlaveProc = 0;

% Set times
DEParams.maxiter  = nsamples/DEParams.NP;
DEParams.maxtime  = 30; % in seconds
DEParams.maxclock = [];
DEParams.VTR  = 1e-7;

% Start differential evolution
[bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
	DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams); %#ok

xout = [bestFctParams.parameter1(1),bestFctParams.parameter2(1),bestFctParams.parameter3(1),bestFctParams.parameter4(1)];
end

function yx = fxval(xs,d,example_idx)
   xs = [xs.parameter1(1),xs.parameter2(1),xs.parameter3(1),xs.parameter4(1)];
   if example_idx == 1
       yx1 = 0; yx2 = 1;
       for jd = 1:d
           yx1 = yx1+(xs(:,jd)+0.1).^2./4;
           yx2 = yx2.*cos((xs(:,jd)+0.1)./sqrt(jd));
        end
       yx = 1+yx1-yx2;
   elseif example_idx == 2  
        yx = 0;
        for jd = 1:d-1
            wj = 1+(xs(:,jd)-1)./4;
            yx = yx + (wj-1).^2.*(1+10.*sin(pi.*wj+1).^2);
        end
        wd = 1+(xs(:,d)-1)./4; w1 = 1+(xs(:,1)-1)./4;
        yx = yx + sin(pi.*w1).^2 + (wd-1).^2.*(1+sin(2*pi.*wd).^2); yx = (yx./1);
   elseif example_idx == 3
        yx = 0;
       for jd = 1:d
           yx = yx + (xs(:,jd)+0.1).^2-10.*cos(2*pi.*(xs(:,jd)+0.1));
       end
       yx = yx + 10.*d; yx = yx./200;
   elseif example_idx == 4
        yx1 = 0; yx2 = 0;  
       for jd = 1:d
           yx1 = yx1+ (xs(:,jd)+0.1).^2;
           yx2 = yx2 + 0.5.*(xs(:,jd)+0.1).*jd;
       end
       yx = yx1 + yx2.^2 + yx2.^4;  yx = yx./100;
    elseif example_idx == 5
   yx1 = 0; yx2 = 0;  
   for jd = 1:d
       yx1 = yx1+ (xs(:,jd)+0.1).^2;
       yx2 = yx2 + cos(2*pi.*(xs(:,jd)+0.1));
   end
   yx = -20.*exp(-0.2.*sqrt(1./d.*yx1)) -exp(1./d.*yx2) + 20 + exp(1);
   yx = yx./100;
   end
end