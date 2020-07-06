function [fx,minfx,qmin] = V(q,x)
    fx = chart(q,x);
    psi0 = chart(q,fz(zeros(2,1)));
    [minfx,qaux] = min(0.5*[norm(chart(-1,x)-psi0)^2 norm(chart(1,x)-psi0)^2]);
    qmin = qaux*2-3;
    fx = 0.5*norm(fx-psi0)^2;
end