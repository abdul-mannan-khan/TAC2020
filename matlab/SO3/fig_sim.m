addpath('../')
k = 1; %controller gain
%% Isomorphism from lie algebra of SO(3) and (R3,x)
D2sk = mf([3 1 27 3]', @(x) zeros(27,3));
Dsk  = mf([3 1 9 3]', @(x) -[S(ee(1,3));S(ee(2,3));S(ee(3,3))],D2sk);
sk   = mf([3 1 3 3]', @(x) S(x),Dsk);
%and its inverse
D2ski= mf([3 3 27 9]', @(x) zeros(27,9));
Dski = mf([3 3 3 9]', @(x) (0.5*eye(3))*[0 0 0 0 0 1 0 -1 0;0 0 -1 0 0 0 1 0 0;0 1 0 -1 0 0 0 0 0],D2ski);
ski  = mf([3 3 3 1]', @(x) [x(3,2)-x(2,3);x(1,3)-x(3,1);x(2,1)-x(1,2)]*0.5,Dski);
%% Basis of the tangent space to S4
%PT(\xs),\PT(\xr)
D2R = mf([3 3 81 9]', @(x) zeros(81,9));
DR  = mf([3 3 9 9]', @(x) eye(9),D2R);
R   = mf([3 3 3 3]', @(x) x,DR);
PT = -kron(eye(3),R)*[S(ee(1,3));S(ee(2,3));S(ee(3,3))];
%PT =kron(eye(3),X(1,1))+sk(X(2:4,1));
%PT = sk(X(2:4,1));
%% Cayley map inverse
D2R = mf([3 3 81 9]', @(x) zeros(81,9));
DR  = mf([3 3 9 9]', @(x) eye(9),D2R);
R   = mf([3 3 3 3]', @(x) x,DR);
Ci = inv(eye(3)+R)*(eye(3)-R);
%% charts
D2R = mf([3 3 81 9]', @(x) zeros(81,9));
DR  = mf([3 3 9 9]', @(x) eye(9),D2R);
R   = mf([3 3 3 3]', @(x) x,DR);
charts{1} = ski(Ci(R));
Dcharts{1} = D(charts{1});
for q = 2:4
    Rshift = axangle(ee(q-1,3),pi/2);
    charts{q} = ski(Ci(Rshift*R));
    Dcharts{q} = D(charts{q});
end
%% Kappa
%\kappa(\theta,\xs,\ws)
DX = mf([21 1 21 21]',@(x) eye(21));
X  = mf([21 1 21 1]',@(x) x,DX);
Kappa = -PT(reshape(X(10:18,1),[3 3]))'*X(1:9,1)-X(19:21,1); 
%% Theta
%\theta(\xr,\xs)
tic
DX = mf([18 1 18 18]',@(x) eye(18));
X  = mf([18 1 18 1]',@(x) x,DX);
for q = 1:4
    TT{q} = Dcharts{q}(reshape(X(10:18,1),[3 3]))'*...
        (charts{q}(reshape(X(10:18,1),[3 3]))-charts{q}(reshape(X(1:9,1),[3 3])));
    DTT{q} = D(TT{q});
end
toc
%% Lyapunov function
%\W{q}(\xr,\xs,\ws,\theta)
tic
DX = mf([30 1 30 30]',@(x) eye(30));
X  = mf([30 1 30 1]',@(x) x,DX);
for q = 1:4
    W{q} = 0.5*(charts{q}(reshape(X(10:18,1),[3 3]))-charts{q}(reshape(X(1:9,1),[3 3])))'*(charts{q}(reshape(X(10:18,1),[3 3]))-charts{q}(reshape(X(1:9,1),[3 3])))+...
        0.5*(X(19:21,1)'*X(19:21,1))+0.5*((TT{q}(X(1:18,1))-X(22:30,1))'*(TT{q}(X(1:18,1))-X(22:30,1)));
    DW{q} = D(W{q});
end
toc
%%
R0 = diag([1 -1 -1]);
y = vec(eye(3));
[t,j,xi] = run(50,vec(R0),[0,0,0]',Kappa,TT,DTT,PT,W,DW,charts,Dcharts);
%%
layout = [1;2]*ones(1,3);
h = create_axis(layout,15);
axes(h(1))
z = 3-xi(:,1:9)*y;
plotarc(t,j,z)
grid on
set(gca,'xticklabel','','ylim',enlarge([0 max(z)],1.1))
ylbl = ylabel('$3-\xs\projt(t)\tp\xr{}$');
ypos = get(ylbl,'position');
set(ylbl,'position',ypos-[2 0 0])
ypos = get(ylbl,'position');
axes(h(2))
z = sqrt(sum(xi(:,10:12).^2,2));
plotarc(t,j,z)
set(gca,'ylim',enlarge([0 max(z)],1.1))
grid on
xlbl = xlabel('$t$ [s]');
set(xlbl,'position',get(xlbl,'position')+[0 -0.05 0]);
ylbl = ylabel('$\norm{\ws\projt(t)}$');
set(ylbl,'position',get(ylbl,'position')*diag([0 1 1])+[ypos(1) 0 0]);

%{
idx = ~logical([diff(j);0]);
N = numel(t);
Wk = zeros(N,1);
Wk2 = zeros(N,1);
uc = zeros(N,1);
ud = zeros(N,1);
gapc = zeros(N,1);
gapd = zeros(N,1);
ttdist = zeros(N,1);
for I = 1:N
    if idx(I)
        xr = vec(eye(3));
        xs = xi(I,1:9)';
        ws = xi(I,10:12)';
        q  = xi(I,13);
        htt= xi(I,14:22)';
        tt = TT{q}([xr;xs]);
        gapc(I) = mu(xs,ws,q,htt,W);
        uc(I) = norm(Kappa([htt;xs;ws]));
        ttdist(I) = norm(htt-tt);
        Wk(I) = W{q}([xr;xs;ws;htt]);
    end
end
plot(t(idx),Wk(idx))

layout = (1:4)';
h = create_axis(layout,15);
axes(h(1))
plot(t,sqrt(sum((repmat(xr',[numel(t),1])-xi(:,1:9)).^2,2)))
grid on
%legend({'continuous control','discontinuous control'})
ylabel('$\norm{\xs-\xr}$')
axes(h(2))
plot(t,uc)
grid on
%legend({'continuous control','discontinuous control'})
ylabel('$\norm{u\projt(t)}$')
axes(h(3))
plot(t,gapc)
%legend({'Continuous control)','Discontinuous control'})
ylabel('$\gap(\xi\projt(t))$')
grid on
axes(h(4))
plot(t,ttdist)
grid on
ylabel('$\norm{\htt\projt(t)-\tts(\xi\projt(t))}')
xlabel('t [s]')
%}

rmpath('../')