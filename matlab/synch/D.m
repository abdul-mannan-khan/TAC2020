function v = D(xi)
    global delta
    x = xi(1:4);
    q = xi(9);
    Vq = charts(x,q)^2;
    Vnq = charts(x,-q)^2;
    if isnan(Vnq) 
        gap = 0;
    else
        gap = Vq-Vnq;
    end
    if isnan(Vq)
        v = 1;
    elseif gap < delta
        v = 0;
    else
        v = 1;
    end