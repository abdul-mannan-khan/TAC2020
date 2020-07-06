addpath('../')
k = 1; %controller gain
%% Isomorphism from lie algebra of SO(3) and (R3,x)
D2sk = mf([3 1 27 3]', @(x) zeros(27,3));
Dsk  = mf([3 1 9 3]', @(x) -[S(ee(1,3));S(ee(2,3));S(ee(3,3))],D2sk);
sk   = mf([3 1 3 3]', @(x) S(x),Dsk);
%% Basis of the tangent space to S4
%PT(\xs),\PT(\xr)
D2X = mf([4 1 16 4]', @(x) zeros(16,4));
DX  = mf([4 1 4 4]', @(x) eye(4),D2X);
X   = mf([4 1 4 1]', @(x) x,DX);
PT = (0.5*eye(4))*([-X(2:4,1)';kron(eye(3),X(1,1))+sk(X(2:4,1))]);
%PT =kron(eye(3),X(1,1))+sk(X(2:4,1));
%PT = sk(X(2:4,1));

%% Stereographic projection
%\charts{q}(\xs),\charts{q}(\xr)
tic
D2X = mf([4 1 16 4]', @(x) zeros(16,4));
DX  = mf([4 1 4 4]', @(x) eye(4),D2X);
X   = mf([4 1 4 1]', @(x) x,DX);
charts{1} = X(2:4,1)*inv(1+X(1,1));
Dcharts{1}= D(charts{1});
charts{2} = X(2:4,1)*inv(1-X(1,1));
Dcharts{2}= D(charts{2});
toc
%% Kappa
%\kappa(\theta,\xs,\ws)
DX = mf([11 1 11 11]',@(x) eye(11));
X  = mf([11 1 11 1]',@(x) x,DX);
Kappa = -PT(X(5:8,1))'*X(1:4,1)-X(9:11,1); 
%% Theta
%\theta(\xr,\xs)
tic
DX = mf([8 1 8 8]',@(x) eye(8));
X  = mf([8 1 8 1]',@(x) x,DX);
TT{1} = Dcharts{1}(X(5:8,1))'*(charts{1}(X(5:8,1))-charts{1}(X(1:4,1)));
TT{2} = Dcharts{2}(X(5:8,1))'*(charts{2}(X(5:8,1))-charts{2}(X(1:4,1)));
DTT{1} = D(TT{1});
DTT{2} = D(TT{2});
toc
%% Lyapunov function
%\W{q}(\xr,\xs,\ws,\theta)
tic
DX = mf([15 1 15 15]',@(x) eye(15));
X  = mf([15 1 15 1]',@(x) x,DX);
W{1} = 0.5*((charts{1}(X(5:8,1))-charts{1}(X(1:4,1)))'*(charts{1}(X(5:8,1))-charts{1}(X(1:4,1))))+...
    0.5*(X(9:11,1)'*X(9:11,1))+0.5*((TT{1}(X(1:8,1))-X(12:15,1))'*(TT{1}(X(1:8,1))-X(12:15,1)));
W{2} = 0.5*((charts{2}(X(5:8,1))-charts{2}(X(1:4,1)))'*(charts{2}(X(5:8,1))-charts{2}(X(1:4,1))))+...
    0.5*(X(9:11,1)'*X(9:11,1))+0.5*((TT{2}(X(1:8,1))-X(12:15,1))'*(TT{2}(X(1:8,1))-X(12:15,1)));
DW{1} = D(W{1});
DW{2} = D(W{2});
toc
%%
layout = ones(1,3);
ha = create_axis(layout,15);
h = [];
for J =-1:2:1
    [t,j,xi] = run(100,[0;1;0;0],[0,0,0]',J,Kappa,TT,DTT,PT,W,DW,charts,Dcharts);
    h = [h plot(t,2*acos(xi(:,5)),'Linewidth',(J+3)/2)];
    hold all
end
set(gca,'ylim',enlarge([0 2*pi],1.1),...
    'ytick',0:pi/2:2*pi,'yticklabel',{'0','\frac{\pi}{2}','\pi','\frac{3\pi}{2}','2\pi'})
grid on
xlabel('$t$ [s]')
ylabel('$2\arccos(\xs_1)$')
[~,~,~,a] = legend({'$\qs=-1$','$\qs=1$'},'position',[6.5e-01   3.7500e-01   0.3   3.6444e-01]);
hold off
rmpath('../')
%{
    %
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
        if 1%idx(I)
            xr = xi(I,1:4)';
            xs = xi(I,5:8)';
            ws = xi(I,9:11)';
            q  = xi(I,12);
            htt= xi(I,13:16)';
            h = (q+1)/2+1;
    %         h2 = (xi(I,32)+1)/2+1;
            tt = TT{h}([xr;xs]);
    %         tt2 = TT{3-h2}([xi(I,1:6) xi(I,26:31) xi(I,23:25)]');
    %         gapd(I) = max([W{h2}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt']')-...
    %             W{3-h2}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt2']'),0]);
            gapc(I) = max([W{h}([xr;xs;ws;htt])-...
                W{3-h}([-xr;xs;ws;htt]),0]);
            uc(I) = norm(Kappa([htt;xs;ws]));
    %         ud(I) = norm(Kappa([tt' xi(I,29:31) xi(I,26:28)]'));
            ttdist(I) = norm(htt-tt);
            Wk(I) = W{h}([xr;xs;ws;htt]);
            Wk2(I) = W{3-h}([-xr;xs;ws;htt]);
        end
    end
%plot(t,[Wk xi(:,end)])
layout = (1:4)';
h = create_axis(layout,15);
axes(h(1))
plot(t,sqrt(sum((xi(:,1:4)-xi(:,5:8)).^2,2)))
%plot(t,[Wk Wk2])
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