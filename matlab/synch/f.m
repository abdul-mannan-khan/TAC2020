function dxi = f(xi)
    x = xi(1:4);
    w = xi(5:8);
    
    u  = law(xi);
    
    dx = PTSn(x)*u;
    dw = u;
    dxi = [dx;dw;0];
end
function fx = PTSn(x)
    fx = eye(4)-blkdiag(x(1:2)*x(1:2)',x(3:4)*x(3:4)');
end