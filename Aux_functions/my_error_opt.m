function y = my_error_opt(x1,x2)
    n = max(size(x1,1),size(x1,2));
    y = norm(reshape(x1,[n,1])-reshape(x2,[n,1]),inf);
end