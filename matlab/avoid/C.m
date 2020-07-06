function out = C(xi)
    global delta
    x = xi(1:3);
    q = xi(4);
    [Vx,minVx] = V(q,x);
    if Vx-minVx <= delta
        out = 1;
    else
        out = 0;
    end
end