addpath('../')
%%
create_functions
%%
[t,j,xi] = run(5,[-1;0;0]/sqrt(1+0^2),[0,0,10]',FX,Kappa,W,TT,DTT,PT,wb,DW,charts,Dcharts,0);
%%
idx = ~logical([diff(j);0]);
N = numel(t);
Wk = zeros(N,1);
uc = zeros(N,1);
ud = zeros(N,1);
gapc = zeros(N,1);
gapd = zeros(N,1);
ttdist = zeros(N,1);
h = zeros(N,1);
h2 = zeros(N,1);
for I = 1:N
    if 1%idx(I)
        h(I) = (xi(I,13)+1)/2+1;
        h2(I) = (xi(I,32)+1)/2+1;
        tt = TT{h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25)]');
        tt2 = TT{3-h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25)]');
        gapd(I) = max([W{h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt']')-...
            W{3-h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt2']'),0]);
        gapc(I) = max([W{h(I)}([xi(I,1:12) xi(I,23:25) xi(I,14:22)]')-...
            W{3-h(I)}([xi(I,1:12) xi(I,23:25) xi(I,14:22)]'),0]);
        uc(I) = norm(Kappa([xi(I,14:22) xi(I,10:12) xi(I,7:9)]'));
        ud(I) = norm(Kappa([tt' xi(I,29:31) xi(I,26:28)]'));
        ttdist(I) = norm(xi(I,14:22)-tt');
        %Wk(I) = norm(xi(I,14:22)'-TT{h}([xi(I,1:12)';xi(I,23:25)']));
    end
end
%%
layout = [1 2 3 4];
hax = create_axis(layout,21,...
    'leftmargin',0.01,'rightmargin',0.01,'bottommargin',0.3);
It = arrayfun(@(x) floor(x*numel(t)/4),1:4);
for I = 1:4
    axes(hax(I))
    Ns = 20;
    [X,Y,Z] = sphere(Ns);
    surf(X,Y,Z,'facecolor','none','edgecolor',0.75*ones(1,3))
    hold all
    hl(1) = plot3(xi(1:It(I),7),xi(1:It(I),8),xi(1:It(I),9),'linewidth',2);
    hl(2) = plot3(xi(1:It(I),26),xi(1:It(I),27),xi(1:It(I),28),'linewidth',2);
    hl(3) = plot3(xi(1:It(I),1),xi(1:It(I),2),xi(1:It(I),3),'linewidth',2);
    hold off
    axis equal
    set(gca,'xtick',[-1 1],'ytick',[-1 1],'ztick',[-1 1]);
    xlabel('$x_1$')
    ylabel('$x_2$')
    zlabel('$x_3$')
    text(0,0,1.5,['$t=T_' num2str(I) '$'],...
        'horizontalalignment','center',...
        'verticalalignment','bottom')
end
[~,~,~,a]= legend(hl,{'$\xs$ (continuous control)','$\xr{}$ (discontinuous control)','reference'},...
     'orientation','horizontal','position',[0.1 0.01 0.8 0.11]);
%%
layout = [1;2;3]*ones(1,6);
hax = create_axis(layout,15,'topmargin',0.05);
axes(hax(1))
plotarc(t,j,[sqrt(sum((xi(:,1:3)-xi(:,7:9)).^2,2)) sqrt(sum((xi(:,1:3)-xi(:,26:28)).^2,2))])
ylim = get(hax(1),'ylim');
line((arrayfun(@(x) t(It(x)),(1:3)')*[1 1])', repmat(ylim,[3,1])',...
    'color','k','linestyle','--');
grid on
[~,~,~,a] = legend({'continuous control','discontinuous control'},...
    'position',[0.3788    0.8141    0.5729    0.0871]);
yl = ylabel('$\norm{\xs\projt(t)-\xr{}\projt(t)}$');
yp = get(yl,'position');
set(yl,'position',[-1 yp(2:3)]);
set(gca,'xticklabel','','ylim',enlarge(get(gca,'ylim'),1.1))
arrayfun(@(x) text(t(It(x)),ylim(2)+0.1*(ylim(2)-ylim(1)),['$T_' num2str(x) '$'],...
    'horizontalalignment','center'),1:4)
axes(hax(2))
plotarc(t,j,[uc,ud])
ylim = get(hax(2),'ylim');
line((arrayfun(@(x) t(It(x)),(1:3)')*[1 1])', repmat(ylim,[3,1])',...
    'color','k','linestyle','--');
grid on
yl = ylabel('$\norm{u\projt(t)}$');
yp = get(yl,'position');
set(yl,'position',[-1 yp(2:3)]);
set(gca,'xticklabel','','ylim',enlarge(get(gca,'ylim'),1.1))
axes(hax(3))
plotarc(t,j,ttdist)
set(gca,'ylim',enlarge(get(gca,'ylim'),1.1))
grid on
yl = ylabel('$\norm{\htt\projt(t)-\theta\projt(t)}$');
yp = get(yl,'position');
set(yl,'position',[-1 yp(2:3)]);
xl = xlabel('$t$ [s]');
xp = get(xl,'position');
set(xl,'position',[xp(1) -diff(get(gca,'ylim'))*0.3 xp(3)]);

%%
tf = find(t >= 0.1,1);
layout = ones(1,3);
hax = create_axis(layout,15,'leftmargin',0,'bottommargin',0.2);
axes(hax(1))
plotarc(t(1:tf),j(1:tf),[gapc(1:tf) gapd(1:tf)])
[~,~,~,a] = legend({'$\gap[tt](\htt\projt(t),\qs\projt(t),\hxi\projt(t))$','$\gap(\qs\projt(t),\xi\projt(t))$'},...
    'position',[0.3481    0.6201    0.6182    0.3754]);
grid on
set(gca,'xlim',[0 t(tf)],...
    'ylim',enlarge([0 2],1.1))
xl = xlabel('$t$ [s]');
xp = get(xl,'position');
ylim = get(gca,'ylim');
set(xl,'position',[xp(1),ylim(1) - 0.2* diff(ylim),xp(3)])
set(gca,'xtick',0:0.02:0.1)
%%
rmpath('../')
%
%figure('units','centimeters','position',[0 0 15 15])
%

% plot(t,xi(:,1:3),'--')
% hold all
% plot(t,xi(:,7:9),':')
% plot(t,xi(:,26:28))
% hold on
%plot(t(idx),sqrt(sum((Wk(idx)-xi(idx,end)).^2,2)))
%plot(t(idx),[Wk(idx) xi(idx,end)])
%plot(t,[sqrt(sum((xi(:,1:3)-xi(:,7:9)).^2,2)) sqrt(sum((xi(:,1:3)-xi(:,26:28)).^2,2))])