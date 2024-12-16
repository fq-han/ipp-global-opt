function [yj] = construct_aux_tt_int(d,zz,nz)
    ys = cell(d,1); yj = cell(d,1);
    if size(zz,2) == 1
        zz = repmat(zz,[1,d]);
    end
    parfor jd = 1:d
        ys{jd} = ones(nz,1);
    end
    yone = tt_tensor(ys);
    parfor jd = 1:d
        yjs = yone;
        yjs{jd} = zz(:,jd)';
        yj{jd} = yjs;
    end
end