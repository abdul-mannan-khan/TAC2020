
k = 1; %controller gain
%% Projection onto the plane tangent to x
D2X = mf([3 1 9 3]', @(x) zeros(9,3));
DX  = mf([3 1 3 3]', @(x) eye(3),D2X);
X   = mf([3 1 3 1]', @(x) x,DX);
PT = eye(3)-X*X';
%% Stereographic projection
D3X = mf([3 1 27 3]', @(x) zeros(27,3));
D2X = mf([3 1 9 3]', @(x) zeros(9,3),D3X);
DX  = mf([3 1 3 3]', @(x) eye(3),D2X);
X   = mf([3 1 3 1]', @(x) x,DX);
C1 = [eye(2),zeros(2,1)]*X;
C2 = inv(1-[0 0 1]*X);
C3 = inv(1+[0 0 1]*X);
charts{1} = C1*C3;
Dcharts{1}= D(charts{1});
charts{2} = C1*C2;
Dcharts{2}= D(charts{2});
%% Right inverse
D2X= mf([2 3 36 6]', @(x) zeros(36,6));
DX = mf([2 3 6 6]', @(x) eye(6), D2X);
X  = mf([2 3 2 3]', @(x) x,DX);
RInv = X'*inv(X*X');
%% Wback
D2X= mf([9 1 81 9]',@(x) zeros(81,9));
DX = mf([9 1 9 9]',@(x) eye(9),D2X);
X  = mf([9 1 9 1]',@(x) x,DX);
% wb{1} = Dcharts{1}(X(7:9,1));
wb{1} = RInv(Dcharts{1}(X(7:9,1))*PT(X(7:9,1)))*Dcharts{1}(X(1:3,1))*PT(X(1:3,1))*X(4:6,1);
wb{2} = RInv(Dcharts{2}(X(7:9,1))*PT(X(7:9,1)))*Dcharts{2}(X(1:3,1))*PT(X(1:3,1))*X(4:6,1);
Dwb{1} = D(wb{1});
Dwb{2} = D(wb{2});
%% Fx
DX = mf([15 1 15 15]',@(x) eye(15));
X  = mf([15 1 15 1]',@(x) x,DX);
FX = [PT(X(1:3,1))*X(4:6,1);X(13:15,1);PT(X(7:9,1))*X(10:12,1)];
%% Kappa
DX = mf([15 1 15 15]',@(x) eye(15));
X  = mf([15 1 15 1]',@(x) x,DX);
Kappa = X(7:9,1)...
    -k*PT(X(13:15,1))'*X(1:3,1)...
    -k*(X(10:12,1)-X(4:6,1));
%% Theta
DX = mf([15 1 15 15]',@(x) eye(15));
X  = mf([15 1 15 1]',@(x) x,DX);
TT{1} = [Dcharts{1}(X(7:9,1))'*(charts{1}(X(7:9,1))-charts{1}(X(1:3,1)));
         wb{1}(X(1:9,1));Dwb{1}(X(1:9,1))*FX(X)];
TT{2} = [Dcharts{2}(X(7:9,1))'*(charts{2}(X(7:9,1))-charts{2}(X(1:3,1)));
         wb{2}(X(1:9,1));Dwb{2}(X(1:9,1))*FX(X)];
DTT{1} = D(TT{1});
DTT{2} = D(TT{2});
%% Lyapunov function
DX = mf([24 1 24 24]',@(x) eye(24));
X  = mf([24 1 24 1]',@(x) x,DX);
Wpsi{1} = 0.5*k*(charts{1}(X(7:9,1))-charts{1}(X(1:3,1)))'*(charts{1}(X(7:9,1))-charts{1}(X(1:3,1)));
Wpsi{2} = 0.5*k*(charts{2}(X(7:9,1))-charts{2}(X(1:3,1)))'*(charts{2}(X(7:9,1))-charts{2}(X(1:3,1)));
Ww{1} = 0.5*(X(10:12,1)-wb{1}(X(1:9,1)))'*(X(10:12,1)-wb{1}(X(1:9,1)));
Ww{2} = 0.5*(X(10:12,1)-wb{2}(X(1:9,1)))'*(X(10:12,1)-wb{2}(X(1:9,1)));
Wtt{1}= 0.5*(X(16:24,1)-TT{1}(X(1:15,1)))'*(X(16:24,1)-TT{1}(X(1:15,1)));
Wtt{2}= 0.5*(X(16:24,1)-TT{2}(X(1:15,1)))'*(X(16:24,1)-TT{2}(X(1:15,1)));
W{1} = Wpsi{1}+Ww{1}+Wtt{1};
W{2} = Wpsi{2}+Ww{2}+Wtt{2};
% W{1} = 0.5*(X(16:24,1)-TT{1}(X(1:15,1)))'*(X(16:24,1)-TT{1}(X(1:15,1)));
% W{2} = 0.5*(X(16:24,1)-TT{2}(X(1:15,1)))'*(X(16:24,1)-TT{2}(X(1:15,1)));
DW{1} = D(W{1});
DW{2} = D(W{2});