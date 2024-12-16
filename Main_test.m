clear all; close all;
%% Test routine for optimiation of benchmark functions
% term_tol, the stopping criteria for |x_k-x^*|, 
% d dimension, delta coefficients for Gibbs measure
% example index, see choose_example, contians various benchamark functions 
rng(10); %random number generator
addpath(genpath('TT_toolbox'), genpath('Optimization_Algorithms'), 'Aux_functions');
term_tol = 1e-2;   examples = [1]; nexamples = 1;  
d =4; Max_eval = 1e6; Maxit = floor(Max_eval/(40*d));  
delta = 0.1; epsf = 1e-1; % initial delta should not be too small which may cause high computational cost in the beginning
xinit = (rand(1,d)-0.5).*6;
 
% for parameters in TT-IPP and MC-IPP, see the paper
%% Settting TT-IPP_Params
TTIPP_Params.term_tol = term_tol; TTIPP_Params.Maxit = 200; TTIPP_Params.d = d; 
TTIPP_Params.epsf = epsf; TTIPP_Params.delta = delta;  
TTIPP_Params.m = 4; TTIPP_Params.Ch = 1e3; TTIPP_Params.refine_gamma = 1.1; %control whether to recompute TT approximation
TTIPP_Params.etap = 2; TTIPP_Params.etam = 0.5;  TTIPP_Params.tk1 = 1;
TTIPP_Params.alphak = 1; TTIPP_Params.tau = min(TTIPP_Params.tk1/2,1/2);  TTIPP_Params.T = 20; 
TTIPP_Params.warmstart = 1; 

%% Setting MC-IPP Params
MCIPP_Params.term_tol = term_tol; MCIPP_Params.Maxit = Maxit; MCIPP_Params.d = d; 
MCIPP_Params.delta = delta; MCIPP_Params.NMC_int = 40*d; NMC_int = 40*d;
MCIPP_Params.max_evals = Max_eval; 
MCIPP_Params.tk2 = 1;
MCIPP_Params.etap = 2; MCIPP_Params.etam = 0.90; MCIPP_Params.rej_rate = 0.8; %rejection_rate
MCIPP_Params.tau = min(MCIPP_Params.tk2/2,1/2);  MCIPP_Params.T = 20; MCIPP_Params.m = 4;
MCIPP_Params.c_delta = 0.9; MCIPP_Params.C_Nint = 1.1; 
MCIPP_Params.alphak = 0.2; MCIPP_Params.alpha_min = 0.2; MCIPP_Params.alpha_max = 0.3;
MCIPP_Params.warmstart = 1; 

for jexample = 1:nexamples  
    example_idx = examples(jexample); [fval,~,xex] = choose_example(1,delta,d,example_idx);
    % [xk1,xk1_hist,errxk1,nsamples1,k1] = TT_IPP_optimization(example_idx,xinit,TTIPP_Params);
    % [xk2,xk2_hist,errxk2,nsamples2,k2] = MC_HJ_prox_optimization(example_idx,delta,d,floor(Max_eval/NMC_int),NMC_int,term_tol,xinit);
    % [xk3,xk3_hist,errxk3,nsamples3,k3] = my_CBO_optimization(example_idx,d, Maxit,term_tol,xinit);
    [xk4,xk4_hist,errxk4,nsamples4,k4] = my_DE_optimization(example_idx,d, Maxit,term_tol,xinit);
    % [xk5,xk5_hist,errxk5,nsamples5,k5] = my_PSO_optimization(example_idx,d, Maxit,term_tol,xinit);
    % [xk7,xk7_hist,errxk7,nsamples7,k7] = MC_IPP_prox(example_idx,xinit,MCIPP_Params);
    % [xk6,xk6_hist,errxk6,nsamples6,k6] = my_simulanneal(example_idx,d,Maxit*(20*d),term_tol,xinit);
    % [xk8,xk8_hist,errxk8,nsamples8,k8] = my_random_search(example_idx,d,Maxit,term_tol,xinit);
end
 
fprintf(['TT-IPP: ', num2str(errxk1(k1)), '; Function evalutaions: ', num2str(nsamples1), '\n']);
fprintf(['MC_IPP: ',  num2str(my_error_opt(xk7,xex)),'; Function evalutaions: ', num2str(nsamples7), '\n']);
fprintf(['HJ-Prox: ',  num2str(my_error_opt(xk2,xex)),'; Function evalutaions: ', num2str(nsamples2), '\n']);
fprintf(['CBO: ',  num2str(my_error_opt(xk3,xex)),'; Function evalutaions: ', num2str(nsamples3), '\n']);
fprintf(['DE: ', num2str(errxk4(k4)),'; Function evalutaions: ', num2str(nsamples4), '\n']);
fprintf(['PSO: ', num2str(errxk5(k5)),'; Function evalutaions: ', num2str(nsamples5), '\n']);
fprintf(['SA: ', num2str(errxk6(k6)),'; Function evalutaions: ', num2str(nsamples6), '\n']);
fprintf(['PS: ',num2str(my_error_opt(xk8,xex)), '; Function evalutaions: ', num2str(nsamples8), '\n']);

%xki is the final optimization value, nsamples is the number of samples
%used in MC type methods, ki is the number of iterations
%errxki is the distance between xk and x^* 
 
