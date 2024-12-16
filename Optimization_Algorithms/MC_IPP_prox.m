function [xk2,xk_hist,errxk2,nsamples,k] = MC_IPP_prox(example_idx,xinit,MCIPP_Params)
    d= MCIPP_Params.d;
    delta = MCIPP_Params.delta;  
    [fval,~,xex] = choose_example(1,delta,d,example_idx);
    [xk2,xk_hist,errxk2,nsamples,k] = MC_IPP_prox_fun(fval,xinit,xex,MCIPP_Params);
end