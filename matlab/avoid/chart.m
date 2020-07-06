function [fx,dfx] = psi(q,x)
    if x(3)*q ~= -1
        fx = [x(1);x(2)/(1+q*x(3))];
        dfx = [1/(1+x(3)*q),-x(2)*q/(1+q*x(3))^2];
        dfx = blkdiag(1,dfx);
    else
        fx = Inf(2,1);
        dfx = NaN(2,3);
    end
end