function fx = Kmn(m,n)
%
% Kmn - Commutation matrix
%
% fx = Kmn(m,n)
%  m,n - are natural numbers
%
% DESCRIPTION
%  Commutation matrices are such that $K_{mn}\vec(A)=\vec(A\tp)$ for each
%  $m\x n$ matrix (see [Section 9.2.1]{HandbookOfMatrices})
%
% PARENTS
%  \vec
% CHILD
%  \newcommand{\Kmn}[2]{K_{#1#2}}

fx = zeros(m*n);
for I=1:m
    for J = 1:n
        H = zeros(m,n);
        H(I,J) = 1;
        fx = fx+kron(H,H');
    end
end