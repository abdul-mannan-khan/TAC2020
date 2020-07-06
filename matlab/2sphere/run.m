function [t,j,xi] = run(T,xs,ws,FX,kappa,W,TT,DTT,PT,wb,DW,charts,Dcharts,noise)    
    xr = [1;0;0];
    wr = [0;1;0];
    dwr= [-1;0;0];
    htt = zeros(9,1);
    q  = 1;
    h = (q+1)/2+1;
    options = odeset('reltol',1e-6);
    TSPAN = [0 T];
    JSPAN = [0 10];
    rule = 1;
    d = 1;
    if xs(3) == q
        W0 = Inf;
    elseif xs(3) == -q
        W0 = 0;
    else
        tt = TT{h}([xr;wr;xs;ws;dwr]);
        W0 = W{h}([xr;wr;xs;ws;dwr;tt]);
    end
    [t,j,xi] = HyEQsolver(@(xi) F(xi,kappa,FX,TT,DTT,wb,PT,DW,charts,Dcharts,noise),...
            @(xi) G(xi,d,W,TT),...
            @(xi) C(xi,d,W,TT),...
            @(xi) D(xi,d,W,TT),...
            [xr;wr;xs;ws;q;htt;dwr;xs;ws;q;W0;xs;ws],TSPAN,JSPAN,rule,options);
end

function v = C(xi,d,W,TT)
    xr = xi(1:3);
    wr = xi(4:6);
    xs = xi(7:9);
    q  = xi(13);
    htt= xi(14:22);
    dwr = xi(23:25);
    if xs(3) == q
        gap = Inf;
    elseif xs(3) == -q
        gap = 0;
    else
        h  = (q+1)/2+1;
        gap = W{h}([xi(1:12);dwr;htt])-W{3-h}([xi(1:12);dwr;htt]);
    end
    xs2= xi(26:28);
    ws2= xi(29:31);
    q2 = xi(32);
    if xs2(3) == q2
        gap2 = Inf;
    elseif xs2(3) == -q2
        gap2 = 0;
    else
        h2  = (q2+1)/2+1;
        tt2 = TT{h2}([xr;wr;xs2;ws2;dwr]);
        tt2p = TT{3-h2}([xr;wr;xs2;ws2;dwr]);
        gap2 = W{h2}([xr;wr;xs2;ws2;dwr;tt2])-W{3-h2}([xr;wr;xs2;ws2;dwr;tt2p]);
    end
    if gap <= d && gap2 <= d
        v = 1;
    else
        v = 0;
    end
end

function xip = G(xi,d,W,TT)
    xip = xi;
    xr = xi(1:3);
    wr = xi(4:6);
    xs = xi(7:9);
    q  = xi(13);
    htt= xi(14:22);
    dwr = xi(23:25);
    if xs(3) == q
        gap = Inf;
    elseif xs(3) == -q
        gap = 0;
    else
        h  = (q+1)/2+1;
        gap = W{h}([xi(1:12);dwr;htt])-W{3-h}([xi(1:12);dwr;htt]);
    end
    xs2= xi(26:28);
    ws2= xi(29:31);
    q2 = xi(32);
    if xs2(3) == q2
        gap2 = Inf;
    elseif xs2(3) == -q2
        gap2 = 0;
    else
        h2  = (xip(32)+1)/2+1;
        tt2 = TT{h2}([xr;wr;xs2;ws2;dwr]);
        tt2p = TT{3-h2}([xr;wr;xs2;ws2;dwr]);
        gap2 = W{h2}([xr;wr;xs2;ws2;dwr;tt2])-W{3-h2}([xr;wr;xs2;ws2;dwr;tt2p]);
    end
    if gap >= d
        xip(13) = -xip(13);
    end
    if gap2 >= d
        xip(32) = -xip(32);
        h2  = (xip(32)+1)/2+1;
        tt2 = TT{h2}([xr;wr;xs2;ws2;dwr]);
        xip(33) = W{h2}([xip(1:6);xip(26:31);xip(23:25);tt2]);
    end
end

function v = D(xi,d,W,TT)
        xr = xi(1:3);
    wr = xi(4:6);
    xs = xi(7:9);
    q  = xi(13);
    htt= xi(14:22);
    dwr = xi(23:25);
    if xs(3) == q
        gap = Inf;
    elseif xs(3) == -q
        gap = 0;
    else
        h  = (q+1)/2+1;
        gap = W{h}([xi(1:12);dwr;htt])-W{3-h}([xi(1:12);dwr;htt]);
    end
    xs2= xi(26:28);
    ws2= xi(29:31);
    q2 = xi(32);
    if xs2(3) == q2
        gap2 = Inf;
    elseif xs2(3) == -q2
        gap2 = 0;
    else
        h2  = (q2+1)/2+1;
        tt2 = TT{h2}([xr;wr;xs2;ws2;dwr]);
        tt2p = TT{3-h2}([xr;wr;xs2;ws2;dwr]);
        gap2 = W{h2}([xr;wr;xs2;ws2;dwr;tt2])-W{3-h2}([xr;wr;xs2;ws2;dwr;tt2p]);
    end
    if gap < d && gap2 < d
        v = 0;
    else
        v = 1;
    end
end
    
    
    