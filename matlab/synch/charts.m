function [fx,dfx] = charts(x,q)
    S = [0 1;-1 0];
    x1 = x(1:2);
    x2 = x(3:4);
    [fx,df1] = stereo([x1'*x2;q*x1'*S*x2]);
    dfx = df1*[x2', x1';q*x2'*S',q*x1'*S];
end
function [fx,dfx] = stereo(x)
    fx = x(1)/(1+x(2));
    dfx = [1/(1+x(2)),-x(1)/(1+x(2))^2];
end