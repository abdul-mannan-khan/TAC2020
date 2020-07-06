function out = vec(X)
%
% out = vec(X)
%  X - mxn matrix 
%  out - mn-dimensional column vector that stacks the columns of X on top of
%   each other
% 
% PARENTS
%  \newcommand{\ceq}{:=}
%  \newcommand{\x}{times}
%  \newcommand{\bmtx}[1]{\begin{bmatrix}#1\end{bmatrix}}
%  \newcommand{\ee}[2]{e_{#1}^{#2}}
%  \newcommand{\R}[1]{\mathbb{R}^{#1}}
%   
%  CHILD
%   \DeclareMathOperator{\vec}{vec}
%   
%  DESCRIPTION
%   The operator $\vec$ is given by 
%   \begin{equation}\label{eq:vec}
%   \vec(M)\ceq\bmtx{M \ee[n]{1}\\ \vdots \\ M \ee[n]{n}}
%   \end{equation}
%   for each $M\in\R{m\x n}$.
%
    out = X(:);