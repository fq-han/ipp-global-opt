function  [xk7,xk7_hist,errxk7,nsamples7,k7] = my_random_search(example_idx,d,Maxit,term_tol,xinit)
    [fval,~,xex] = choose_example(1,1,d,example_idx);
    nsamples7 = Maxit*(20*d);
    % option.ConstraintTolerance = 1e-20; option.MeshTolerance = 1e-20; option.StepTolerance = 1e-20;
    ffval = zeros(nsamples7,1);
    zj = (rand(nsamples7,d)-0.5).*10;
    for jsample = 1:nsamples7
    ffval(jsample) = fval(zj(jsample,:));
    end
    [~,idx] = min(ffval);
    xk7 = zj(idx,:);
    k7 = 1;
    % nsamples7 = output.funccount; k7 = output.iterations;
    errxk7 = (my_error_opt(xk7,xex)); xk7_hist = xk7;
end