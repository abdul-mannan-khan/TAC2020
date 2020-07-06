function [out,d] = fx(x)
    global z0 ee
    out = z0+x(2:3)*(exp(x(1))+ee);
    d = [x(2:3)*exp(x(1)) (exp(x(1))+ee)*eye(2)];
end