function [t,j,xi] = run(xi0)
    TSPAN = [0 10];
    JSPAN = [0 10];
    options = odeset('reltol',1e-6);
    rule  = 1;
    [t,j,xi] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
end


function out = D(xi)
    global delta
    x = xi(1:3);
    q = xi(4);
    [Vx,minVx] = V(q,x);
    if Vx-minVx < delta
        out = 0;
    else
        out = 1;
    end
end
function out = g(xi)
    out = xi;
    out(4) = -out(4);
end



    