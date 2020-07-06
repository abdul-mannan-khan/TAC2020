function [fv,v] = mu(xs,ws,q,htt,W)
    xr = eye(3);
    qs = [];
    for I = 1:4
        if I == 1
            Rshift = eye(3);
        else
            Rshift = axangle(ee(I-1,3),pi/2);
        end
        if det(Rshift*reshape(xs,[3 3])+eye(3)) ~= 0
            qs = [qs I];
        end
    end
    
    if all(qs ~= q)
        fv = Inf;
        [~,v]  = min(arrayfun(@(h) W{h}([vec(xr);xs;ws;htt]),qs));
        v = qs(v);
    else
        fv = 0;
        v  = q;
        for p = qs 
            aux = W{q}([vec(xr);xs;ws;htt])-W{p}([vec(xr);xs;ws;htt]);
            if aux > fv
                fv = aux;
                v  = p;
            end
        end
    end
end