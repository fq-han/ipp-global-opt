function ttfw = add_weight_tt(ttf,hz,d)
    core_ttf = core(ttf);
    if length(hz) == 1
        hz = hz.*ones(1,d);
    end
    parfor jd = 1:d
        core_ttf{jd} = core_ttf{jd}.*hz(jd);
    end
    ttfw = tt_tensor(core_ttf);  
end