%
create_functions
%
x0 = [-1;0;0]/sqrt(1+0^2);
w0 = [0,0,10]';
noise = 0;
[t,j,xi] = run(5,x0,w0,FX,Kappa,W,TT,DTT,PT,wb,DW,charts,Dcharts,noise);
%
idx = ~logical([diff(j);0]);
N = numel(t);
Wc = zeros(N,3);
Wd = zeros(N,2);
h = zeros(N,1);
h2 = zeros(N,1);
for I = 1:N
    if 1%idx(I)
        h(I) = (xi(I,13)+1)/2+1;
        h2(I) = (xi(I,32)+1)/2+1;
        tt = TT{h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25)]');
        tt2 = TT{3-h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25)]');
        Wc(I,1) = Wpsi{h(I)}([xi(I,1:12) xi(I,23:25) xi(I,14:22)]');
        Wd(I,1) = Wpsi{h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt2']');
        Wc(I,2) = Ww{h(I)}([xi(I,1:12) xi(I,23:25) xi(I,14:22)]');
        Wd(I,2) = Ww{h2(I)}([xi(I,1:6) xi(I,26:31) xi(I,23:25) tt2']');
        Wc(I,3) = Wtt{h(I)}([xi(I,1:12) xi(I,23:25) xi(I,14:22)]');
    end
end
%
%%
layout = [ones(4) 2*ones(4);
          9*ones(1,8);
          3*ones(2,6) 4*ones(2,2);
          5*ones(2,6) 6*ones(2,2);
          7*ones(2,6) 8*ones(2,2)];
hax = create_axis(layout,21,...
    'leftmargin',0.05,'rightmargin',0.05,'bottommargin',0.05,'topmargin',0.05,...
    'innerymargin',0.01,'innerxmargin',0.01);
set(hax(9),'visible','off')
It = arrayfun(@(x) floor(x*numel(t)/3),1:3);
It = It([1 end]);
for I = 2:-1:1
    axes(hax(I))
    Ns = 20;
    [X,Y,Z] = sphere(Ns);
    surf(X,Y,Z,'facecolor','none','edgecolor',0.75*ones(1,3))
    hold all
    hl(1) = plot3(xi(1:20:It(I),7),xi(1:20:It(I),8),xi(1:20:It(I),9),'linewidth',2,'linestyle',':');
    hl(2) = plot3(xi(1:20:It(I),26),xi(1:20:It(I),27),xi(1:20:It(I),28),'linewidth',2);
    hl(3) = plot3(xi(1:20:It(I),1),xi(1:20:It(I),2),xi(1:20:It(I),3),'linewidth',2,'linestyle','--');
    hold off
    axis equal
    set(gca,'xtick',[-1 1],'ytick',[-1 1],'ztick',[-1 1]);
    xlabel('$x_1$')
    ylabel('$x_2$')
    text(0,0,1.5,['$t=T_' num2str(I) '$'],...
        'horizontalalignment','center',...
        'verticalalignment','bottom')
end
zlabel('$x_3$')
[~,~,~,a]= legend(hl,{'Section~\ref{sub:continuous}','Section~\ref{sub:hcontrol}','$\xr{\qs}(t)$'},...
     'orientation','horizontal','position',[0.05   0.95    0.9    0.05]);
%
%layout = (1:3)'*ones(1,4);
%hax = create_axis(layout,10);
axes(hax(3))
plotarc(t,j,Wc(:,1),[],[],{':','linewidth',2})
hold all
colororder = get(gca,'colororder');
haux = plotarc(t,j,Wd(:,1),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6});
hold off
grid on
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'xticklabel','','ylim',ylim)
text(t(end)-0.01*diff(t([1 end])),ylim(2)-0.01*diff(ylim),...
    '$\displaystyle\frac{1}{2}\norm{\psi{\qs}(\xs)-\psi{\qs}(\xr{\qs})}^2$',...
    'interpreter','none','horizontalalignment','right','verticalalignment','top')

%[~,~,~,a] = legend({'continuous control','discontinuous control'},...
%   'position',[0.4403    0.6311    0.2927    0.0711]);
line((arrayfun(@(x) t(It(x)),(1:2)')*[1 1])', repmat(ylim,[2,1])',...
    'color','k','linestyle','--');
arrayfun(@(x) text(t(It(x)),ylim(2)+0.1*(ylim(2)-ylim(1)),['$T_' num2str(x) '$'],...
    'horizontalalignment','center'),1:2)
axes(hax(4))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),1),[],[],{':','linewidth',2})
hold on
plotarc(t(1:It(1)),j(1:It(1)),Wd(1:It(1),1),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
set(gca,'xticklabel','','yaxislocation','right','xlim',[0 t(It(1))])
grid on
%
axes(hax(5))
plotarc(t,j,Wc(:,2),[],[],{':','linewidth',2})
hold all
plotarc(t,j,Wd(:,2),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
grid on
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'xticklabel','','ylim',ylim)
line((arrayfun(@(x) t(It(x)),(1:2)')*[1 1])', repmat(ylim,[2,1])',...
    'color','k','linestyle','--');
text(t(end)-0.01*diff(t([1 end])),ylim(2)-0.01*diff(ylim),...
    '$\displaystyle\frac{1}{2}\norm{\ws-\Tmc[\qs](\xs,\xr{\qs})\wr{\qs}}^2$',...
    'interpreter','none','horizontalalignment','right','verticalalignment','top')

axes(hax(6))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),2),[],[],{':','linewidth',2})
hold on
plotarc(t(1:It(1)),j(1:It(1)),Wd(1:It(1),2),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
set(gca,'xticklabel','','yaxislocation','right','xlim',[0 t(It(1))])
grid on
%
axes(hax(7))
plotarc(t,j,Wc(:,3),[],[],{':','linewidth',2})
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'ylim',ylim)
text(t(end)-0.01*diff(t([1 end])),ylim(2)-0.01*diff(ylim),...
    '$\displaystyle\frac{1}{2}\norm{\htt-\theta(\zs)}^2$',...
    'interpreter','none','horizontalalignment','right','verticalalignment','top')
xlbl = xlabel('$t$','interpreter','none');
set(xlbl,'units','normalized')
set(xlbl,'position',get(xlbl,'position')-[0 .1 0])
grid on
line((arrayfun(@(x) t(It(x)),(1:2)')*[1 1])', repmat(ylim,[2,1])',...
    'color','k','linestyle','--');
axes(hax(8))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),3),[],[],{':','linewidth',2})
set(gca,'yaxislocation','right','xlim',[0 t(It(1))])
grid on
xlbl = xlabel('$t$','interpreter','none');
set(xlbl,'units','normalized')
set(xlbl,'position',get(xlbl,'position')-[0 0.1 0])
%{
layout = [ones(4) 2*ones(4) 3*ones(4) 4*ones(4);
          11*ones(1,16);
          5*ones(3,12) 6*ones(3,4);
          7*ones(3,12) 8*ones(3,4);
          9*ones(3,12) 10*ones(3,4)];
hax = create_axis(layout,21,...
    'leftmargin',0.05,'rightmargin',0.05,'bottommargin',0.05,'topmargin',0.05,...
    'innerymargin',0.01,'innerxmargin',0.01);
set(hax(11),'visible','off')
It = arrayfun(@(x) floor(x*numel(t)/4),1:4);
for I = 1:4
    axes(hax(I))
    Ns = 20;
    [X,Y,Z] = sphere(Ns);
    surf(X,Y,Z,'facecolor','none','edgecolor',0.75*ones(1,3))
    hold all
    hl(1) = plot3(xi(1:20:It(I),7),xi(1:20:It(I),8),xi(1:20:It(I),9),'linewidth',2,'linestyle',':');
    hl(2) = plot3(xi(1:20:It(I),26),xi(1:20:It(I),27),xi(1:20:It(I),28),'linewidth',2);
    hl(3) = plot3(xi(1:20:It(I),1),xi(1:20:It(I),2),xi(1:20:It(I),3),'linewidth',2,'linestyle','--');
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
[~,~,~,a]= legend(hl,{'Controller of Section~\ref{sub:continuous}','Controller of Section~\ref{sub:hcontrol}','reference'},...
     'orientation','horizontal','position',[0.05   0.95    0.9    0.05]);
%
%layout = (1:3)'*ones(1,4);
%hax = create_axis(layout,10);
axes(hax(5))
plotarc(t,j,Wc(:,1),[],[],{':','linewidth',2})
hold all
colororder = get(gca,'colororder');
haux = plotarc(t,j,Wd(:,1),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6});
hold off
grid on

ylabel('$\displaystyle\frac{1}{2}\norm{\psi{\qs}(\xs)-\psi{\qs}(\xr{\qs})}^2$','interpreter','none')
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'xticklabel','','ylim',ylim)
%[~,~,~,a] = legend({'continuous control','discontinuous control'},...
%   'position',[0.4403    0.6311    0.2927    0.0711]);
line((arrayfun(@(x) t(It(x)),(1:4)')*[1 1])', repmat(ylim,[4,1])',...
    'color','k','linestyle','--');
arrayfun(@(x) text(t(It(x)),ylim(2)+0.1*(ylim(2)-ylim(1)),['$T_' num2str(x) '$'],...
    'horizontalalignment','center'),1:4)
axes(hax(6))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),1),[],[],{':','linewidth',2})
hold on
plotarc(t(1:It(1)),j(1:It(1)),Wd(1:It(1),1),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
set(gca,'xticklabel','','yaxislocation','right','xlim',[0 t(It(1))])
grid on
%
axes(hax(7))
plotarc(t,j,Wc(:,2),[],[],{':','linewidth',2})
hold all
plotarc(t,j,Wd(:,2),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
grid on
ylabel('$\displaystyle\frac{1}{2}\norm{\ws-\wback{\qs}}^2$','interpreter','none')
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'xticklabel','','ylim',ylim)
line((arrayfun(@(x) t(It(x)),(1:4)')*[1 1])', repmat(ylim,[4,1])',...
    'color','k','linestyle','--');
axes(hax(8))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),2),[],[],{':','linewidth',2})
hold on
plotarc(t(1:It(1)),j(1:It(1)),Wd(1:It(1),2),[],[],{'color',colororder(2,:),'linewidth',2},...
    {'color',colororder(2,:),'linestyle','--','marker','*','markerfacecolor',colororder(2,:),'markersize',6})
hold off
set(gca,'xticklabel','','yaxislocation','right','xlim',[0 t(It(1))])
grid on
%
axes(hax(9))
plotarc(t,j,Wc(:,3),[],[],{':','linewidth',2})
ylabel('$\displaystyle\frac{1}{2}\norm{\htt-\theta(\zs)}^2$','interpreter','none')
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'ylim',ylim)
xlabel('$t$','interpreter','none')
grid on
line((arrayfun(@(x) t(It(x)),(1:4)')*[1 1])', repmat(ylim,[4,1])',...
    'color','k','linestyle','--');
axes(hax(10))
plotarc(t(1:It(1)),j(1:It(1)),Wc(1:It(1),3),[],[],{':','linewidth',2})
set(gca,'yaxislocation','right','xlim',[0 t(It(1))])
grid on
xlabel('$t$','interpreter','none')
%}