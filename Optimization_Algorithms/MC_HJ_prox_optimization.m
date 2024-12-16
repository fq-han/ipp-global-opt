function [xk2,xk_hist,errxk2,nsamples,k] = MC_HJ_prox_optimization(example_idx,delta,d,Maxit,NMC_int,term_tol,xinit)
    tk2 = 1;  xk_hist = zeros(Maxit,d);
    etap = 2; etam = 0.1; warmstart =1;
    alphak = 1-1/2*sqrt(etam);  
    tau = min(tk2/2,1/2);  T = 4*tk2; 
    f2 = 100;     errxk2 = zeros(Maxit,1);
    [fval,finval,xex] = choose_example(1,delta,d,example_idx);
    tol_qk = term_tol*1e-2;                
    
    if_shrink_delta = 0;
    nx = 1; qk2 = ones(nx,d);  
    xk2 = reshape(xinit,[1,d]);  finit = fval(xk2);   fmin = 100;  bestsample = xinit;
    if warmstart == 1
    Mcvd = 0; Mcgvd = 0;  
    for j = 1:NMC_int
       zj = (rand(1,d)-0.5).*6; 
       fzj = fval(zj);
       if fzj<fmin
            bestsample = zj; fmin = fzj;
        end
        Mcvd = Mcvd + exp(-fzj./delta);
        Mcgvd = Mcgvd + zj.*exp(-fzj./delta);
    end
    Mcvd = Mcvd./NMC_int;
    Mcgvd = Mcgvd./NMC_int;
    xk2 = Mcgvd./(Mcvd*ones(1,d));
    else
    xk2 = xinit;
    end
    errxk2(1) = my_error_opt(xk2,xex);

    for k = 2:Maxit
        %% By MC integration
        if norm(qk2,inf)>tol_qk
            Mcvd = 0; Mcgvd = 0;  
            for j = 1:NMC_int
                zj = sqrt(delta*tk2).*randn(1,d) + xk2; 
                fzj = fval(zj);
                if fzj<fmin
                    bestsample = zj; fmin = fzj;
                end
                Mcvd = Mcvd + exp(-fzj./delta);
                Mcgvd = Mcgvd + zj.*exp(-fzj./delta);
            end
            Mcvd = Mcvd./NMC_int;
            Mcgvd = Mcgvd./NMC_int;
            proxtf2 = Mcgvd./(Mcvd*ones(1,d));
            qk2p = qk2;
            qk2 = (xk2-proxtf2);    
            xk2 = xk2 - alphak*qk2;
            tk2 = timestep(tk2,qk2,qk2p,tau,T,etam,etap);
        else
            fprintf('HJ-MAD trapped');
            xk2 = bestsample;
            nsamples = k*NMC_int;
            return;
        end
        xk_hist(k,:) = xk2;
        f2k = fval(xk2);
        if if_shrink_delta == 1
            if f2k > f2 - min(1,finit)*1e-3/k && delta > 1e-7
               delta = delta/2;  
               fprintf(['MC shrink delta, delta = ', num2str(delta),'\n']);
            end
        end
        f2 = f2k;
        errxk2(k) = my_error_opt(xk2,xex);
        if errxk2(k) < term_tol
           nsamples = k*NMC_int;
            xk2 = bestsample;
           return;
        end
        if isnan(errxk2(k))
            nsamples = k*NMC_int;
            xk2 = bestsample;
            return;
        end
        if fval(xk2) < fmin
           bestsample = xk2;
           fmin = fval(xk2);
        end
    end
    xk2 = bestsample;
    nsamples = k*NMC_int;
end