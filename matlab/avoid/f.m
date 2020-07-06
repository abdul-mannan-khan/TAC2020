function dxi = f(xi)
    x = xi(1:3);
    q = xi(4);
    w = xi(5:6);
    [fpsi,dpsi] = chart(q,x);
    psi0 = chart(q,fz(zeros(2,1)));
    [~,dfz] = fz(fx(x));
    u = -dfz'*dpsi'*(fpsi-psi0)-w;
    dx = dfz*w;
    dw = u;
    dxi = [dx;0;dw];
end