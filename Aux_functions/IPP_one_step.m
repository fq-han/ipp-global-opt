function [xk,tk,qk] =  IPP_one_step(xk,zz,ttf,qk,delta,tk,yj,tau,T,alphak,d,nx,etam,etap)
    proxtf1 = zeros(nx,d);   Kdts = cell(d,1);
    zz = reshape(zz,[length(zz),1]);
    tic;
    for jx = 1:nx
        for jd = 1:d
            Kdts{jd}  = exp(-(zz-xk(jx,jd)).^2./(delta*tk*2));
        end
        Kdt = tt_tensor(Kdts); 
        vdelta = dot(Kdt,ttf);  
        tic;
        parfor jd = 1:d
            gvdeltajd = dot(Kdt,times(yj{jd},ttf));
            proxtf1(jx,jd) = (gvdeltajd)./(vdelta);
        end
        toc;
    end
    qkp = qk;
    qk = (xk-proxtf1);
    xk = xk - alphak*qk;
    tk  = timestep(tk,qk,qkp,tau,T,etam,etap);
end

 