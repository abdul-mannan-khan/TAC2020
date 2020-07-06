global z0 ee delta 
z0 = [1;0];
ee = 0.5;
delta = 1;

z_init = [2;0];
x0 = fz(z_init);
z = [];
w = [];
for I = 1:2
    xi0 = [x0;I*2-3;zeros(2,1)];
    [t,j,xi] = run(xi0);
    for J = 1:numel(t)
        z(:,J,I) = fx(xi(J,1:3)');
        w(:,J,I) = xi(J,5:6)';
    end
end
%%
Nr = 100;
Ntt= 101;
rlim = [-3 log(4)];
xlim = [-0.5 2.5];
h = create_axis(ones(5,4),15);
set(gca,'xlim',xlim)
axis equal
hold all
hc = draw_circle('center',[1;0],'radius',0.5,'N',40);
set(hc,'color','k')
step = ceil(numel(t)/10);
colors = get(gcf,'defaultaxescolororder');
linestyles = {'-','--'};
s = 0.8;
for I = 1:2
    hc(I+1) = fig_fset(I*2-3,Nr,Ntt,rlim,colors(I,:)*(1-s)+ones(1,3)*s,linestyles{I},1);
end
hc(4) = fig_fset(0,Nr,Ntt,rlim,colors(3,:)*(1-s)+ones(1,3)*s,linestyles{1},2);
for I = 1:2
    hc(I+4) = plot(z(1,:,I)',z(2,:,I)','color',colors(I,:));
    quiver(z(1,1:step:end,I),z(2,1:step:end,I),w(1,1:step:end,I),w(2,1:step:end,I),'color',colors(I,:));
end
set(hc(5),'linewidth',2)
grid on
xlabel('$z_1$')
ylabel('$z_2$')
[~,~,~,a] = legend(hc,{'Obstacle','$C_{-1}\minus C_1$','$C_{1}\minus C_{-1}$','$C_{-1}\cap C_{1}$',...
    '$\qs(0,0)=-1$','$\qs(0,0)=1$'},...
    'position',[1.5169e-01   7.3519e-01   3.0515e-01   2.4412e-01]);
quiver(z0(1),z0(2),ee,0,'color','k')
text(z0(1)+ee/2,0,'$\epsilon$','horizontalalignment','c','verticalalignment','bottom')
hold off
set(h,'ylim',get(h,'ylim')+0.3)