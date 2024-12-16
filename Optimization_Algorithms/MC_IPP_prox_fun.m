function [xk2,xk_hist,errxk2,nsamples,k] = MC_IPP_prox_fun(fval,xinit,xex,MCIPP_Params)
    tk2 = MCIPP_Params.tk2;  
    max_evals = MCIPP_Params.max_evals; d= MCIPP_Params.d;
    delta = MCIPP_Params.delta; Maxit = MCIPP_Params.Maxit;
    etap = MCIPP_Params.etap; etam = MCIPP_Params.etam; rej_rate = MCIPP_Params.rej_rate; %rejection_rate
    tau = MCIPP_Params.tau;  T = MCIPP_Params.T; m = MCIPP_Params.m;
    c_delta = MCIPP_Params.c_delta; C_Nint = MCIPP_Params.C_Nint;
    alphak = MCIPP_Params.alphak; alpha_min = MCIPP_Params.alpha_min; alpha_max = MCIPP_Params.alpha_max;
    NMC_int = MCIPP_Params.NMC_int;  term_tol = MCIPP_Params.term_tol;
 
    errxk2 = zeros(Maxit,4); 
    jhist = 1;
    bestsample = fval(xinit); fmin = fval(xinit);
    if MCIPP_Params.warmstart == 1
    Mcvd = 0; Mcgvd = 0;
    for j = 1:NMC_int
        zj = (rand(1,d)-0.5).*6;
        fzj = fval(zj);
        if fzj<fmin
            bestsample = zj;
            fmin = fzj;
        end
        expfval = exp(-fzj/delta);
        Mcvd = Mcvd + expfval;
        Mcgvd = Mcgvd + zj.*expfval;
    end
    xk2 = Mcgvd./Mcvd;
    else
    xk2 = xinit;
    end
    xk_hist = zeros(Maxit,d);
    nx = 1; qk2 = ones(nx,d);  
    mean_error_decay = -1; nsamples = 0; finit = fval(xk2);
    xk2p = xk2;
    for k = 1:Maxit
        [prox_pt,bestsample,fmin] = proximal_point(fval,delta,tk2,nx,d,xk2,1:d,NMC_int,bestsample,fmin);
        if nsamples > max_evals
            xk2 = bestsample;
            fprintf('MCIPP hits maximum number of samples');
            k = k-1;
            return;
        end
       
        xk2 = alphak*prox_pt+(1-alphak)*xk2p;
        nsamples = nsamples + NMC_int;
        if nnz(isnan(xk2))>0
           xk2 = xk2p;
           tk2 = tk2*etam; 
           xk_hist(k,:) = xk2;
           errxk2(k,1) = my_error_opt(xk2,xex);   errxk2(k,2) = fval(xk2);
           errxk2(k,3) = delta; errxk2(k,4) = NMC_int; errxk2(k,5) = tk2;
           continue;
        end
        if fval(xk2) > fval(xk2p) && rand(1) < rej_rate
           xk2 = xk2p;
        else
            xk2p = xk2; xk_hist(jhist,:) = xk2; jhist = jhist + 1;
            qk2p = qk2;
            qk2 = xk2 - prox_pt;
            tk2 = timestep(tk2,qk2,qk2p,tau,T,etam,etap);
        end
        errxk2(k,1) = my_error_opt(xk2,xex);   errxk2(k,2) = fval(xk2);
        errxk2(k,3) = delta; errxk2(k,4) = NMC_int; errxk2(k,5) = tk2;
        if k > m+1
           mean_error_decay = errxk2(k,2) - max(errxk2(k-m:k-1,2));
        end
        if mean_error_decay>-max(min(1,finit),1e-8)*1e-3/k && delta > 1e-3 
           delta = delta*c_delta;  
           NMC_int =  floor(C_Nint*NMC_int);
           alphak = max(alpha_min,c_delta*alphak);
           fprintf(['MC shrink delta, delta = ', num2str(delta),'\n']);
        else
            alphak = min(alphak/c_delta,alpha_max);
        end
        % my_error_opt(xk2,xex)
        if fval(xk2) < fmin
           bestsample = xk2;
           fmin = fval(xk2);
        end
        if errxk2(k,1) < term_tol
            xk2 = bestsample;
           return;
        end
    end
end

function [prox_pt,bestsample,fmin] = proximal_point(fval,delta,tk2,nx,d,xk2,idset,NMC_int,bestsample,fmin)
    Mcvd = 0; Mcgvd = 0;
    for j = 1:NMC_int         
        zj = sqrt(delta*tk2).*randn(nx,d) + xk2; 
        fzj = fval(zj);
        if fzj<fmin
            bestsample = zj; fmin = fzj;
        end
        expfval = exp(-fzj./delta);
        Mcvd = Mcvd + expfval;
        Mcgvd = Mcgvd + zj(idset).*expfval;
    end
    prox_pt = Mcgvd./Mcvd;
end

