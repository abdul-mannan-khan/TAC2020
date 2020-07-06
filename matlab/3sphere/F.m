function dxi = F(xi,kappa,TT,DTT,PT,DW,charts,Dcharts)
    xr = xi(1:4);
    xs = xi(5:8);
    ws = xi(9:11);
    q  = xi(12);
    htt= xi(13:16);
    h  = (q+1)/2+1;
    tt = TT{h}([xr;xs]);
    %
    dxr = zeros(4,1);
    dxs = PT(xs)*ws;
    dws = kappa([htt;xs;ws]);
    dq = 0;    
    dot_tt = DTT{h}([xr;xs])*[dxr;dxs];
    dhtt = dot_tt-(htt-tt)+dxs;
    %dW = -norm(ws)^2-norm(tt-htt)^2;
    %dW = DW{h}([xr;xs;ws;htt])*[dxr;dxs;dws;dhtt];
    %dW = (charts{h}(xs)-charts{h}(xr))'*Dcharts{h}(xs)*dxs+ws'*dws+(tt-htt)'*(dot_tt-dhtt);
    %dW = tt'*dxs+ws'*dws+(tt-htt)'*(dot_tt-dhtt);
    dW = ws'*PT(xs)'*(tt-htt)+(tt-htt)'*(dot_tt-dhtt)-ws'*ws;
    dxi = [dxr;dxs;dws;dq;dhtt;dW];
end