function [t,j,xi] = run(T,xs,ws,kappa,TT,DTT,PT,W,DW,charts,Dcharts)    
    xr = eye(3);
    htt = zeros(9,1);
    q  = 1;
    options = odeset('reltol',1e-6);
    TSPAN = [0 T];
    JSPAN = [0 10];
    rule = 1;
    d = 1;
    if q == 1
        Rshift = eye(3);
    else
        Rshift = axangle(ee(q-1,3),pi/2);
    end
    if det(Rshift*reshape(xs,[3 3])+eye(3)) == 0
        W0 = Inf;
    else
        W0 = W{q}([vec(xr);xs;ws;htt]);
    end
    [t,j,xi] = HyEQsolver(@(xi) F(xi,kappa,TT,DTT,PT,DW,charts,Dcharts),...
            @(xi) G(xi,d,W),...
            @(xi) C(xi,d,W,TT),...
            @(xi) D(xi,d,W,TT),...
            [xs;ws;q;htt;W0],TSPAN,JSPAN,rule,options);
end

function v = C(xi,d,W,TT)
    xs = xi(1:9);
    ws = xi(10:12);
    q  = xi(13);
    htt= xi(14:22);
    %tt = TT{q}([vec(xr);xs]);
    gap = mu(xs,ws,q,htt,W);
    if gap <= d 
        v = 1;
    else
        v = 0;
    end
end

function xip = G(xi,d,W)
    xr = eye(3);
    xs = xi(1:9);
    ws = xi(10:12);
    q  = xi(13);
    htt= xi(14:22);
    [gap,qp] = mu(xs,ws,q,htt,W);
    if gap >= d
        xip = [xs;ws;qp;htt;W{qp}([vec(xr);xs;ws;htt])];
    else
        xip = [xs;ws;q;htt;W{q}([vec(xr);xs;ws;htt])];
    end
end

function v = D(xi,d,W,TT)
    xs = xi(1:9);
    ws = xi(10:12);
    q  = xi(13);
    htt= xi(14:22);
    %tt = TT{q}([vec(xr);xs]);
    gap = mu(xs,ws,q,htt,W);
    if gap < d 
        v = 0;
    else
        v = 1;
    end
end


                
            