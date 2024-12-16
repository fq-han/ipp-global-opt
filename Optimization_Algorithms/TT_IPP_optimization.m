function [xk1,xk_hist,errxk1,nsamples,k] = TT_IPP_optimization(example_idx,xinit,TTIPP_Params)
    %% TT IPP optimization for benchmark functions
    % create mesh first
    d = TTIPP_Params.d; delta = TTIPP_Params.delta; epsf = TTIPP_Params.epsf;
    Lz =5; nz =51; zz = linspace(-Lz,Lz,nz)'; 
    hz = zz(2)-zz(1);
    zz = repmat(zz,[1,d]); hz = hz.*ones(1,d);
    
    % choose example to optimize
    [fval,finval,xex] = choose_example(zz,delta,d,example_idx);
    [yj] = construct_aux_tt_int(d,zz,nz);         
    % TT cross optimization on a uniform mesh using tt-cross algorithm
    ttf_visitedMap = containers.Map; 
    [ttf,~,~,~,~,nsamples,ttf_visitedMap] = exp_greedy2_cross_interp(nz*ones(1,d), finval, fval, epsf,0,delta,'nswp',10, 'check_dictionary', ttf_visitedMap);

    TTIPP_Params.nsamples = nsamples;  TTIPP_Params.hz = hz; TTIPP_Params.nz = nz; 
    [xk1,xk_hist,errxk1,nsamples,k]= TT_IPP_optimization_fun(ttf,zz,yj,finval,fval,xinit,xex,TTIPP_Params,ttf_visitedMap);
end
