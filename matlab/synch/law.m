function fx = law(xi)
    x = xi(1:4);
    w = xi(5:8);
    q = xi(9);
    [psi,dpsi] = charts(x,q);  
    fx = Pi(x)'*dpsi'*(1-psi);
end
function fx = Pi(x)
    x1 = x(1:2);
    x2 = x(3:4);
    fx = eye(4)-blkdiag(x1*x1',x2*x2');
end