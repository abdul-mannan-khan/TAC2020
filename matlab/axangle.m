function fx = axangle(u,tt)
    fx = eye(3)+sin(tt)*S(u)+(1-cos(tt))*S(u)^2;
end