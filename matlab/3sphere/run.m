function [t,j,xi] = run(T,xs,ws,q,kappa,TT,DTT,PT,W,DW,charts,Dcharts)    
    xr = [-q;0;0;0];
    htt = zeros(4,1);
    h = (q+1)/2+1;
    options = odeset('reltol',1e-6);
    TSPAN = [0 T];
    JSPAN = [0 10];
    rule = 1;
    d = 1;
    if xs(1) == q
        W0 = Inf;
    elseif xs(1) == -q
        W0 = 0;
    else
        W0 = W{h}([xr;xs;ws;htt]);
    end
    [t,j,xi] = HyEQsolver(@(xi) F(xi,kappa,TT,DTT,PT,DW,charts,Dcharts),...
            @(xi) G(xi),...
            @(xi) C(xi,d,W,TT),...
            @(xi) D(xi,d,W,TT),...
            [xr;xs;ws;q;htt;W0],TSPAN,JSPAN,rule,options);
end

function v = C(xi,d,W,TT)
    xr = xi(1:4);
    xs = xi(5:8);
    ws = xi(9:11);
    q  = xi(12);
    h  = (q+1)/2+1;
    htt= xi(13:16);
    tt = TT{h}([xr;xs]);
    if xs(1) == q
        gap = Inf;
    elseif xs(1) == -q
        gap = 0;
    else
        %\W{q}(\xr,\xs,\ws,\theta)
        h  = (q+1)/2+1;
        gap = W{h}([xr;xs;ws;htt])-W{3-h}([-xr;xs;ws;htt]);
    end
    
    if gap <= d 
        v = 1;
    else
        v = 0;
    end
end

function xip = G(xi)
    xip = xi;
    xip(12) = -xip(12);
    xip(1:4)= -xip(1:4);
end

function v = D(xi,d,W,TT)
    xr = xi(1:4);
    xs = xi(5:8);
    ws = xi(9:11);
    q  = xi(12);
    h  = (q+1)/2+1;
    htt= xi(13:16);
    tt = TT{h}([xr;xs]);
    if xs(1) == q
        gap = Inf;
    elseif xs(1) == -q
        gap = 0;
    else
        h  = (q+1)/2+1;
        gap = W{h}([xr;xs;ws;htt])-W{3-h}([-xr;xs;ws;htt]);
    end
    
    if gap < d 
        v = 0;
    else
        v = 1;
    end
end