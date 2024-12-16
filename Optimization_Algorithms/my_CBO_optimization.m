% Consensus based optimization (CBO)
%
% This function performs CBO.
% 
% 
% [vstar_approx] = CBO(E, grad_E, parametersCBO, V0)
% 
% input:    E             = objective function E (as anonymous function)
%           grad_E        = gradient of objective function E (as anonymous function)
%           parametersCBO = suitable parameters for CBO
%                         = [T, dt, N, lambda, gamma, learning_rate, sigma, alpha]
%               - T       = time horizon
%               - dt      = time step size
%               - N       = number of particles
%               - lambda  = consensus drift parameter
%               - gamma   = gradient drift parameter
%               - l._r.   = learning rate associated with gradient drift
%               - sigma   = exploration/noise parameter
%               - alpha   = weight/temperature parameter alpha
%           V0            = initial position of the particles
%           
% output:   vstar_approx  = approximation to vstar
%

function [vstar_approx, xk_hist, errxk, evalcnt, k] = my_CBO_optimization(example_idx,d, Maxit, term_tol,xinit)
    [fun,~,xex] = choose_example(1,1,d,example_idx); 
    N = d*40; V0 = reshape(xinit,[d,1])*ones(1,N) +  (randn(d,N));
     % get parameters
    fmin = 10; bestsample = xinit;
    alpha = 1000;
    dt = 0.01;
    sigma = sqrt(1.6);
    lambda = 1;
    evalcnt = 0;
    
    % initialization
    V = V0;
    
    errxk = zeros(Maxit,1); 
    xk_hist = zeros(Maxit,d);
    for k = 1:Maxit
        
        % % CBO iteration
        % compute current consensus point v_alpha
        [v_alpha,evalcnt,bestsample,fmin] = my_compute_valpha(fun, alpha, V, evalcnt,bestsample,fmin);
    
        % position updates of one iteration of CBO
        V = my_CBO_update(fun, dt,lambda,sigma, v_alpha, V);
        errxk(k) = my_error_opt(v_alpha,xex);
        xk_hist(k,:) = v_alpha;
        if errxk(k) < term_tol
           fprintf('CBO converges');
            vstar_approx = v_alpha;
            return
        end
        if nnz(isnan(v_alpha))>0
           fprintf('CBO has NAN value');
           vstar_approx = v_alpha;
           % vstar_approx = bestsample;
            return
        end
        if fun(v_alpha')<fmin
            bestsample = v_alpha; fmin = fun(v_alpha');
        end
        % if k > 2 && norm(xk_hist(k,:)-xk_hist(k-1,:),inf) < term_tol*1e-2
        %     fprintf('CBO trapped');
        %     vstar_approx = v_alpha;
        %     return;
        % end
    end
    [v_alpha,evalcnt,bestsample,fmin] = my_compute_valpha(fun, alpha, V, evalcnt,bestsample,fmin);
    vstar_approx = v_alpha;
    if fun(v_alpha')>fun(reshape(bestsample,[1,d]))
        vstar_approx = bestsample;
    end
end

function [v_alpha, evalcnt,bestsample,fmin] = my_compute_valpha(fun, alpha, V, evalcnt,bestsample,fmin)
    % energies of the individual particles
    Es = zeros(1,size(V,2));
    for jsamples = 1:size(V,2)
        fVj = fun(V(:,jsamples)');
        Es(jsamples) = fVj;
        if fVj<fmin
            bestsample = V(:,jsamples)'; fmin = fVj;
        end
    end
    
    % minimal energy among the individual particles
    Emin = min(Es);
    
    % computation of current empirical consensus point v_alpha
    w_alpha = exp(-alpha*(Es-Emin));
    v_alpha = sum((V.*w_alpha),2);
    v_alpha = 1/sum(w_alpha)*v_alpha;
    evalcnt = evalcnt + size(V,2);
end