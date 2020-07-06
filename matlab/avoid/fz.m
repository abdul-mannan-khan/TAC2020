function [fx,dfx] = fz(z)
    global z0 ee
    fx = [log(norm(z-z0)-ee);
          (z-z0)/norm(z-z0)];
    dfx= [1/(norm(z-z0)-ee)*(z-z0)'/norm(z-z0);PT((z-z0)/norm(z-z0))/norm(z-z0)];
end

function out = PT(x)
    out = eye(2)-x*x';
end