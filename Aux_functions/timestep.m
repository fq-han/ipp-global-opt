function tk = timestep(tk,qk,qk1,tau,T,etam,etap)
%how to choose the time step t 
%qk1 is the value of qk in the previous iteration
    theta1 = 1/4; theta2= 3/4; epsilon =0.2;
    if norm(qk)<=(theta1*norm(qk1)+epsilon)
        tk = min(etap*tk,T);
    elseif norm(qk)>(theta2*norm(qk1)+epsilon)
        tk = tk;
    else
        tk = max(etam*tk,tau);
    end
end