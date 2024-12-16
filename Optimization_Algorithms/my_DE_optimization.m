function [xk, xk_hist,xkerror, evalcnt, k] = my_DE_optimization(example_idx,d, Maxit, term_tol, xinit)
  %test for differential evolution
    [fun,~,xex] = choose_example(1,1,d,example_idx);
    objFctHandle = @(Params,xs)fxval(xs,d,fun);
    
    % Define parameter names, ranges and quantization:
    
    % 1. column: parameter names
    % 2. column: parameter ranges
    % 3. column: parameter quantizations
    % 4. column: initial values (optional)
    
    Lz = 10;  
    paramDefCell{1,1} = append('parameter');
    paramDefCell{1,2} = ones(d,1)*[-Lz,Lz];
    paramDefCell{1,3} = ones(d,1)*0.01;
    paramDefCell{1,4} = reshape(xinit,[d,1]);
    % Set initial parameter values in struct objFctParams 
    objFctParams.paramDefCell{1,1} = ones(d,1)*4;
    
    % Set single additional function parameter
    objFctSettings = 100;
    
    % Get default DE parameters
    DEParams = getdefaultparams;
    
    % Set number of population members (often 10*D is suggested) 
    DEParams.NP = 40*d;
    
    % Do not use slave processes here. If you want to, set feedSlaveProc to 1 and
    % run startmulticoreslave.m in at least one additional Matlab session.
    DEParams.feedSlaveProc = 0;
    
    % Set times
    DEParams.maxiter  = Maxit;
    % DEParams.maxtime  = 30; % in seconds
    DEParams.maxclock = [];
    DEParams.VTR  = -inf;
    DEParams.displayResults  = 0;
    DEParams.saveHistory = 0;
    DEParams.infoPeriod = 0;
    DEParams.infoIterations = 0;
    DEParams.nochangeiter =200;

    % Start differential evolution
    [bestmem, bestval, bestFctParams, k, evalcnt, xkerror, xk_hist,timeOverFileName] = my_differentialevolution(...
	    DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams,xex,term_tol); %#ok
    xk = bestFctParams.parameter;
    mbdelete(timeOverFileName);
end


function yx = fxval(xs,d, fun)
   xx = xs.parameter';
   % for j = 1:d
   %     % parameterj = append('parameter',num2str(j));
   %     xx(j) = xs.parameterj;
   % end
   yx = fun(xx);
end