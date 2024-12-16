function [xk1,xk_hist,errxk1,nsamples,k]= TT_IPP_optimization_fun(ttf,zz,yj,finval,fval,xinit,xex,TTIPP_Params,ttf_visitedMap)
    %% TT-IPP optimization for function f,0
    % fval is the function handel, finval is the function value on mesh zz,
    % ttf is exp(-f/delta) for initial detla
    % definition of parameters, see paper

    term_tol = TTIPP_Params.term_tol; Maxit = TTIPP_Params.Maxit; d = TTIPP_Params.d; hz = TTIPP_Params.hz;   nz = TTIPP_Params.nz; 
    epsf = TTIPP_Params.epsf; nsamples = TTIPP_Params.nsamples; delta = TTIPP_Params.delta;

    xk_hist = zeros(Maxit,d);  
    tk1 = TTIPP_Params.tk1; 
    etap = TTIPP_Params.etap; etam = TTIPP_Params.etam; m = TTIPP_Params.m;
    alphak = TTIPP_Params.alphak;  tau = TTIPP_Params.tau;  T = TTIPP_Params.T; 
    Ch = TTIPP_Params.Ch; refine_gamma = TTIPP_Params.refine_gamma; %control whether to recompute TT approximation
    errxk1 = zeros(Maxit,2);    
    nx = 1; qk1 = ones(nx,d);
    
    %% using gibbs measure corresponds to f to generate an intial guess
    if TTIPP_Params.warmstart == 1
        gttf = zeros(1,d);
        for jd = 1:d
            gttf(jd) = mytt_dot(tt_ones(nz,d),times(yj{jd},ttf),hz);
        end
        xk1 = (gttf)./mytt_dot(tt_ones(nz,d),ttf,hz);
    else
        xk1 = xinit;
    end
    bestsample = xk1;
    xk_hist(1,:) = xk1;

    finit = fval(xk1); mean_error_decay = -1; fmin = fval(xk1);
    for k = 1:Maxit
        %% By HJ-prox-iteration
        % itk = itk+1;
        [xk1,tk1,qk1] = IPP_one_step(xk1,zz,ttf,qk1,delta,tk1,yj,tau,T,alphak,d,nx,etam,etap);
        f1k = fval(xk1);

        f1 = f1k;
        errxk1(k,1) =  my_error_opt(xk1,xex); 
        errxk1(k,2) = f1; errxk1(k,3) =delta; errxk1(k,4) =epsf; errxk1(k,5) = nz;
        errxk1(k,6) = tk1;
        xk_hist(k+1,:) = xk1;
        if k > m+1
            mean_error_decay = errxk1(k,2) - max(errxk1(k-m:k-1,2));
        end

        if  mean_error_decay>-max(min(1,finit),1e-8)*1e-3/k && delta > 1e-4 && epsf > 1e-7
           delta = delta/2; ttf = times(ttf,ttf); ttf = round(ttf,epsf/1e2, 10); 
           if min(hz) > Ch*delta^refine_gamma
               hz = hz./2; Lz = zz(end)- zz(1); 
               if size(zz,2) == 1
                   zz = repmat(zz,[1,d]); hz = repmat(zz,[1,d]);
               end
               if nz <100*Lz % too many points cause the program crack down
                   nz = (2*nz-1);
               end
               ozz = zz; zz = zeros(nz,d);
               for jd = 1:d
                   [~,xk1_indx] = min(abs(xk1(jd)-ozz(:,jd)));
                   mxk1d = max(-Lz,ozz(xk1_indx,jd)-hz(jd)/2*(nz-1)); 
                   Mxk1d = min(Lz,ozz(xk1_indx,jd)+hz(jd)/2*(nz-1));
                   zz(:,jd) = linspace(mxk1d,Mxk1d,nz)';
               end
               hz = zz(2,:)-zz(1,:); epsf = epsf/2; 
               [ttf,~,~,~,~,nsamples,ttf_visitedMap] = exp_greedy2_cross_interp(nz*ones(1,d), finval, fval, epsf,nsamples,delta,'nswp',10, 'check_dictionary', ttf_visitedMap);
               [yj] = construct_aux_tt_int(d,zz,nz);
           end
        end

         if fval(xk1) < fmin
           bestsample = xk1;
           fmin = fval(xk1);
        end
        if errxk1(k,1) < term_tol
            fprintf('Termination value achieved');
            xk1 = bestsample;
            break;
        end
        if nnz(isnan(errxk1(k,1)))>0
            xk1 = bestsample;
            break;
        end
        if nsamples > 1e6
            fprintf('Max evaluations achieved');
            xk1 = bestsample;
            break;
        end
        % fval(xk1)
    end  
end

function [xk,tk,qk] =  IPP_one_step(xk,zz,ttf,qk,delta,tk,yj,tau,T,alphak,d,nx,etam,etap)
%% one step IPP in tensor format
    proxtf1 = zeros(nx,d);   Kdts = cell(d,1);
    hz = zz(2,:)-zz(1,:);
    for jx = 1:nx
        parfor jd = 1:d
            Kdts{jd}  = exp(-(zz(:,jd)-xk(jx,jd)).^2./(delta*tk*2));
        end
        Kdt = tt_tensor(Kdts); 
        vdelta = mytt_dot(Kdt,ttf,hz);  
        parfor jd = 1:d
            gvdeltajd = mytt_dot(Kdt,times(yj{jd},ttf),hz);  
            proxtf1(jx,jd) = (gvdeltajd)./(vdelta);
        end
    end
    qkp = qk;
    qk = (xk-proxtf1);
    xk = xk - alphak*qk;
    tk  = timestep(tk,qk,qkp,tau,T,etam,etap);
end
