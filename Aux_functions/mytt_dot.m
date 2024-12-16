function [p] = mytt_dot(tt1,tt2,h)
%Dot  product of two TT tensors
%   [PR]=DOT(TT1,TT2) -- dot product of two TT-tensors
%
%   [PR]=DOT(TT1,TT2, DO_QR) if DO_QR==true is specified, perform the 
%   left-to-right QRs of TT1,TT2
%   before the scalar product. It increases the  accuracy in some cases.



d=tt1.d; 
r1=tt1.r; r2=tt2.r; ps1=tt1.ps; ps2=tt2.ps;
n=tt1.n;
core1=tt1.core; core2=tt2.core;

% If the first indices are not ones
p=eye(r1(1)*r2(1));
p = reshape(p, r1(1)*r2(1)*r1(1), r2(1));

for i=1:d
    cr1=core1(ps1(i):ps1(i+1)-1);
    cr2=core2(ps2(i):ps2(i+1)-1);
    cr2=reshape(cr2,[r2(i), n(i)*r2(i+1)]);
    
    p = p*cr2; % size r11*r21*r1-, n*r2+
    p = reshape(p,r1(1)*r2(1), r1(i)*n(i), r2(i+1));
    p = permute(p, [1, 3, 2]);
    p = reshape(p, r1(1)*r2(1)*r2(i+1), r1(i)*n(i));
    
    cr1=reshape(cr1,[r1(i)*n(i), r1(i+1)]);
    
    p = p*conj(cr1).*h(i); % size r11*r12*r2+, r1+
    p = reshape(p, r1(1)*r2(1), r2(i+1), r1(i+1));
    p = permute(p, [1, 3, 2]);
    p = reshape(p, r1(1)*r2(1)*r1(i+1), r2(i+1));
    if i > 4 && min(abs(p(:))) < 1e-100 
        p = p./1e-300;
        % fprintf('rescale p in dot product');
    end
    if i >4 && max(abs(p(:))) > 1e100 
        p = p./1e300;
        % fprintf('rescale p in dot product');
    end
end
end
