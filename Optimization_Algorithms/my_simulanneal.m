function  [xk6,xk6_hist,errxk6,nsamples6,k6] = my_simulanneal(example_idx,d, Maxit,term_tol,xinit)
    [fval,~,xex] = choose_example(1,1,d,example_idx);
    options.Verbosity = 0;  
    [xk6,~,nsamples6,errxk6] = anneal(fval, xinit, Maxit,options,xex,term_tol);
    xk6_hist = xk6; k6 = nsamples6; 
end