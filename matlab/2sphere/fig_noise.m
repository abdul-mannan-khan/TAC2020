create_functions
%
x0 = [-1;0;0]/sqrt(1+0^2);
w0 = [0,0,10]';
noise = 1;
[t,j,xi] = run(14,x0,w0,FX,Kappa,W,TT,DTT,PT,wb,DW,charts,Dcharts,noise);
%
idx = ~logical([diff(j);0]);
N = numel(t);
prtbation = zeros(N,3);
for I = 1:N
    [~,aux] = F(xi(I,:)',Kappa,FX,TT,DTT,wb,PT,DW,charts,Dcharts,noise);
    prtbation(I,:) = aux';
end
%%
layout = (1:2)'*ones(1,3);
hax = create_axis(layout,10);
axes(hax(1));
colors = get(gca,'colororder');
plotarc(t,j,sqrt(sum((xi(:,1:3)-xi(:,7:9)).^2,2)),[],[],{'Color',colors(1,:)});
hold on
plotarc(t,j,sqrt(sum((xi(:,1:3)-xi(:,26:28)).^2,2)),[],[],{'Linewidth',2,'Color',colors(2,:)});
hold on
plotarc(t,j,sqrt(sum((xi(:,1:3)-xi(:,34:36)).^2,2)),[],[],{'--','Color',colors(3,:)});
h1 = findobj(gca,'Type','line','Color',colors(1,:));
h2 = findobj(gca,'Type','line','Color',colors(2,:));
h3 = findobj(gca,'Type','line','Color',colors(3,:));
ylim = enlarge([0 2],1.1);
set(gca,'xticklabel','','ylim',ylim)
ylabel('$\norm{\xs\projt(t)-\xr{}\projt(t)}$','interpreter','none')
grid on
axes(hax(2))
plot(t,sqrt(sum(prtbation.^2,2)))
ylim = enlarge(get(gca,'ylim'),1.1);
set(gca,'ylim',ylim)
xlabel('$t$')
ylabel('$\norm{2\PT(\xs)\grad[\xs]\V{}(\xs,\xr{})}$','interpreter','none')
grid on
axes(hax(1))
[~,~,~,a] = legend([h1(1);h2(1);h3(1)], {'Section~\ref{sub:continuous}',...
    'Section~\ref{sub:hcontrol}',...
    '\cite[Lemma~8]{Bullo1999}'},...
    'position',[.5075    0.4056    0.4585    0.1755]);
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
[~,~,~,a]= legend(hl,{'$\xs$ (continuous control)','$\xr{}$ (discontinuous control)','reference'},...
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