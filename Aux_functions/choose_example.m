function [fval,finval,xex] = choose_example(zz,delta,d,example_idx)
% 1 is Griewank function; 2 is Levy function; 3 is RASTRIGIN FUNCTION
% 4 is ZAKHAROV function 5 is ACKLEY 1 FUNCTION
% 6 Deflected Corrugated Spring function
% 7 Rosenbrock function
% 8 brown function, 9 exponential, 10 Trid function 
% 11 Discus function, 12 Cosine mixture
% 13 Alphine
% 14 biology task DNA optimization, see "Iterative Power Algorithm for GlobalOptimization with Quantics Tensor Trains"
% 15 schaffer 01, 16 Powell (works for dimension up to 16, higher has overflow error), 17 Drop wave
% 18 Bohachevsky 19 Schaffer 02
    
    fval = @(xs)fxval(xs,d,example_idx);
    finval = @(zs)fval(zz(zs));


    if example_idx == 0
       xex = zeros(d,1)-0.1;
    elseif example_idx == 1
       xex = zeros(d,1)-0.1;
    elseif example_idx == 2
        xex = ones(d,1);
   elseif example_idx == 3
        xex = zeros(d,1)-0.1;
   elseif example_idx == 4
        xex = zeros(d,1)-0.1;
   elseif example_idx == 5
        xex = zeros(d,1)-0.1;
    elseif example_idx == 6
        xex = zeros(d,1);
    elseif example_idx == 7
        xex = ones(d,1);
    elseif example_idx == 8
        xex = zeros(d,1);
    elseif example_idx == 9
        xex = zeros(d,1);
   elseif example_idx == 10
        xex = zeros(d,1);
   elseif example_idx == 11
        xex = zeros(d,1);
   elseif example_idx == 12
        xex = zeros(d,1);
    elseif example_idx == 13
        xex = zeros(d,1);
    elseif example_idx == 14
        xex = -1.*ones(d,1);
    elseif example_idx == 15
        xex = zeros(d,1);
    elseif example_idx == 16
        xex = zeros(d,1);
    elseif example_idx == 17
        xex = zeros(d,1);
    elseif example_idx == 18
        xex = zeros(d,1);
    elseif example_idx == 19
        xex = zeros(d,1);
    end
end

function yx = fxval(xs,d,example_idx)
    if example_idx == 0
       yx1 = 0; yx2 = 1;
       for jd = 1:d
           yx1 = yx1+abs(xs(:,jd)+0.1).^(1/2)./4;
           yx2 = yx2.*cos((xs(:,jd)+0.1)./sqrt(jd));
        end
       yx = 1+yx1-yx2;
    elseif example_idx == 1
       yx1 = 0; yx2 = 1;
       for jd = 1:d
           yx1 = yx1+(xs(:,jd)+0.1).^2./4000;
           yx2 = yx2.*cos((xs(:,jd)+0.1)./sqrt(jd));
        end
       yx = (1+yx1-yx2); 
    elseif example_idx == 2
        yx = 0;
        for jd = 1:d-1
            wj = 1+(xs(:,jd)-1)./4;
            wj1 = 1+(xs(:,jd+1)-1)./4;
            yx = yx + (wj-1).^2.*(1+10.*sin(pi.*wj1).^2);
        end
        wd = 1+(xs(:,d)-1)./4; w1 = 1+(xs(:,1)-1)./4;
        yx = yx + sin(pi.*w1).^2 + (wd-1).^2 + yx; yx = (yx);
   elseif example_idx == 3
        yx = 0;
       for jd = 1:d
           yx = yx + (xs(:,jd)+0.1).^2-10.*cos(2*pi.*(xs(:,jd)+0.1));
       end
       yx = yx + 10.*d; yx = yx;
   elseif example_idx == 4
        yx1 = 0; yx2 = 0;  
       for jd = 1:d
           yx1 = yx1+ (xs(:,jd)+0.1).^2;
           yx2 = yx2 + 0.5.*(xs(:,jd)+0.1).*jd;
       end
       yx = yx1 + yx2.^2 + yx2.^4;  yx = yx;
    elseif example_idx == 5
       yx1 = 0; yx2 = 0;  
       for jd = 1:d
           yx1 = yx1+ (xs(:,jd)+0.1).^2;
           yx2 = yx2 + cos(2*pi.*(xs(:,jd)+0.1));
       end
       yx = -20.*exp(-0.2.*sqrt(1./d.*yx1)) -exp(1./d.*yx2) + 20 + exp(1);
       yx = yx;
    elseif example_idx == 6
        yx = 0;
        xsum = sum(xs.^2,2);
        yx = yx + 0.1.*xsum - cos(5*sqrt(xsum));
        yx = yx + 1;
   elseif example_idx == 7
        yx = 0;
        for jd = 1:(d-1)
            term1 = 100 * (xs(:, jd+1) - (xs(:, jd)).^2).^2;
            term2 = (1 - (xs(:, jd))).^2;
            yx = yx + term1 + term2;
        end
        yx = yx;
    elseif example_idx == 8
        yx = 0;
        for jd = 1:d-1
            yx = yx + (xs(:,jd).^2).^(xs(:,jd+1).^2+1)+(xs(:,jd+1).^2).^(xs(:,jd).^2+1);
        end
        yx = yx;
    elseif example_idx == 9
        yx = 0;
        for jd = 1:d
            yx = yx + xs(:,jd).^2;
        end
        yx = -exp(-1/2.*(yx).*d)+1; yx = yx;
    elseif example_idx == 10
        yx = (xs(:,1)-1+d).^2;
        for jd = 2:d
            yx = yx + (xs(:,jd)-1+jd*(d+1-jd)).^2 - (xs(:,jd-1)+(jd-1)*(d+1-jd+1)).*(xs(:,jd)+jd*(d+1-jd));
        end
        yx = yx+d*(d+4)*(d-1)/6; yx = yx;
    elseif example_idx == 11
        yx = 0;
        N = 6000; yx = xs(:,1).^2*N;
        for jd = 2:d
            yx = yx + xs(:,jd).^2;
        end
        yx = yx;
    elseif example_idx == 12
        yx = 0;
        for jd = 1:d
            yx = -0.1*cos(5*pi*xs(:,jd)./20)-(xs(:,jd)./20).^2 + yx;
        end
        yx = yx + 0.1*d; yx = yx./2;
    elseif example_idx == 13
        yx = 0;
        for jd = 1:d
            yx = yx + abs(xs(:,jd).*sin(xs(:,jd))+0.1.*xs(:,jd));
        end
        yx = yx;
    elseif example_idx == 14
        yx = 0;
        for jd = 1:d
            yx = yx + 0.429.*xs(:,jd)-1.126.*xs(:,jd).^2-0.143.*xs(:,jd).^3+0.563.*xs(:,jd).^4;
        end
       yx = yx+0.8490*d;
    elseif example_idx== 15
        yx = 0.5*(d-1);
        for jd = 1:d-1
            sx = xs(:,jd).^2 + xs(:,jd+1).^2;
            yx = yx + (sin(sqrt(sx)).^2-0.5)./(1+0.001.*sx).^2;
        end
        yx = yx;
    elseif example_idx == 16
        yx = 0;
        for jd = 1:floor(d/4)
            yx = ((xs(:,4*jd-3)+10.*xs(:,4*jd-2)).^2 + 5.*(xs(:,4*jd-1)-xs(:,4*jd)).^2 + (xs(:,4*jd-2)-2.*xs(:,4*jd-1)).^4 + 10.*(xs(:,4*jd-3)-xs(:,4*jd)).^4)+yx;
        end
        yx = yx./10;
    elseif example_idx == 17
        yx1 = 0; % Accumulator for the squared sum
        % Compute the sum of squares of each point
        for jd = 1:d
            yx1 = yx1 + xs(:, jd).^2;
        end
        % Calculate the Drop Wave function
        yx = -(1+cos(12.*sqrt(yx1)))./ (2 + 0.5 * yx1);
        yx = yx + 1;
    elseif example_idx == 18
        yx = 0; 
        for jd = 1:d-1
            yx = yx + xs(:, jd).^2 + 2*xs(:,jd+1).^2 - 0.3.*cos(3*pi*xs(:,jd)) - 0.4.*cos(4*pi.*xs(:,jd+1))+0.7;
        end
    elseif example_idx == 19
        yx = 0.5*(d-1);
        for jd = 1:d-1
            dx = xs(:,jd).^2 - xs(:,jd+1).^2;
            sx = xs(:,jd).^2 + xs(:,jd+1).^2;
            yx = yx + (sin((dx)).^2-0.5)./(1+0.001.*sx).^2;
        end
        yx = yx;
    end
end

