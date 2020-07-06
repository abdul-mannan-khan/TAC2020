function dxi = F(xi,kappa,TT,DTT,PT,DW,charts,Dcharts)
    xr = eye(3);
    xs = xi(1:9);
    ws = xi(10:12);
    q  = xi(13);
    htt= xi(14:22);
    tt = TT{q}([vec(xr);xs]);
    %
    dxr = zeros(9,1);
    dxs = PT(reshape(xs,[3 3]))*ws;
    dws = kappa([htt;xs;ws]);
    dq = 0;    
    dot_tt = DTT{q}([vec(xr);xs])*[dxr;dxs];
    dhtt = dot_tt-(htt-tt)+dxs;
    %dW = -norm(ws)^2-norm(tt-htt)^2;
    %dW = DW{h}([xr;xs;ws;htt])*[dxr;dxs;dws;dhtt];
    %dW = (charts{h}(xs)-charts{h}(xr))'*Dcharts{h}(xs)*dxs+ws'*dws+(tt-htt)'*(dot_tt-dhtt);
    %dW = tt'*dxs+ws'*dws+(tt-htt)'*(dot_tt-dhtt);
    dW = ws'*PT(reshape(xs,[3 3]))'*(tt-htt)+(tt-htt)'*(dot_tt-dhtt)-ws'*ws;
    dxi = [dxs;dws;dq;dhtt;dW];
end