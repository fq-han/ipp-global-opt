function [Vq] = tt_interp_nd(zz,nx,d,f_tt,xs)
    Vq = zeros(nx,1); zidx = zeros(d,1);  
    if d == 3
        for jx = 1:nx
            xj = xs(jx,:); 
            for jd = 1:d
                diffz = abs(xj(jd)-zz); [~,zidx(jd)] = min(diffz);
            end
            Vq(jx) = f_tt(zidx(1),zidx(2),zidx(3));
        end
    elseif d == 6
        for jx = 1:nx
            xj = xs(jx,:); 
            for jd = 1:d
                diffz = abs(xj(jd)-zz); [~,zidx(jd)] = min(diffz);
            end
            Vq(jx) = f_tt(zidx(1),zidx(2),zidx(3),zidx(4),zidx(5),zidx(6));
        end
    else
        for jx = 1:nx
            xj = xs(jx,:); 
            for jd = 1:d
                diffz = abs(xj(jd)-zz); [~,zidx(jd)] = min(diffz);
            end
            tmp = num2cell(zidx);
            Vq(jx) = f_tt(tmp{:});
        end
    end
    % Vq = zeros(nx,1); Vj = zeros(2,2,2,2,2,2); 
    % for j = 1:nx
    %     x1j = xs(j,1); x2j = xs(j,2); x3j = xs(j,3);
    %     x4j = xs(j,4); x5j = xs(j,5); x6j = xs(j,6);
    %     zj = value_look([x1j,x2j,x3j,x4j,x5j,x6j],zz,nz);
    %     zj(zj>nz) = nz; zj(zj<2) = 2;
    %     zj = [zj;zj-1];
    %     [izj1,izj2,izj3,izj4,izj5,izj6] = ndgrid(zj(:,1),zj(:,2),zj(:,3),zj(:,4),zj(:,5),zj(:,6)); 
    %     zj1 = zz(izj1);  zj2 = zz(izj2);  zj3 = zz(izj3); 
    %     zj4 = zz(izj4);  zj5 = zz(izj5);  zj6 = zz(izj6); 
    % 
    %     dist = (zz(izj1)-x1j).^2+(zz(izj2)-x2j).^2+(zz(izj3)-x3j).^2+(zz(izj4)-x4j).^2+(zz(izj5)-x5j).^2+(zz(izj6)-x6j).^2;
    %     [~,minI] = min(dist(:));
    %     Vq(j) = f_tt(izj1(minI),izj2(minI),izj3(minI),izj4(minI),izj5(minI),izj6(minI));
    %     % Vq(j) = interpn(zj1,zj2,zj3,zj4,zj5,zj6,Vj,x1j,x2j,x3j,x4j,x5j,x6j,'spline');
    %  end
end

% function Vq  = my_interpolation(zj1,zj2,zj3,Vj,x1,x2,x3)
%     dist = sqrt((zj1-x1).^2+(zj2-x2).^2+(zj3-x3).^2);
%     Vq = sum(Vj(:)./dist(:))./sum(1./dist(:));
% end
% 
% function zj = value_look(xj,zz,nz)
%     halfnz = (nz+1)/2; zj = zeros(size(xj));
%     for j = 1:length(xj)
%         if xj(j) > 0
%            ij = 1; zj(j) = halfnz;
%            while ij < halfnz
%                  if xj(j)>zz(ij)
%                     zj(j) = ij;
%                     break;
%                  end
%                  ij = ij + 1;
%            end
%         elseif xj(j) < 0
%            ij = halfnz; zj(j) = nz;
%            while ij < nz
%                  if xj(j)>zz(ij)
%                     zj(j) = ij;
%                     break;
%                  end
%                  ij = ij + 1;
%            end
%         end
%     end
% end