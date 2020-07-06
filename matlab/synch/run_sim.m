global delta
delta  = 2;
TSPAN = [0 10];
JSPAN = [0 10];
rule  = 1;

%{
x0 = [];

for I = 1:2
    aux = rand(2,1)*2-1;
    aux = aux/norm(aux);
    x0 = [x0;aux];
end
x0 = [x0;rand(4,1)*2-1;1];
%}
x0 = [1;0;0;-1;zeros(4,1);1];

[t,j,x] = HyEQsolver(@f,@g,@C,@D,x0,TSPAN,JSPAN,rule);
%%
h = create_axis(1,11,'leftmargin',0,'bottommargin',0,'rightmargin',0);
hp= plot([x(:,1) x(:,3)],[x(:,2) x(:,4)],'linewidth',2);
%axis equal
tt = linspace(0,2*pi,100);
line(cos(tt),sin(tt),'color','k','linestyle','--')
set(gca,'xlim',[-1.15 1.15],'ylim',[-1.15 1.15],'box','on','xtick','',...
    'ytick','','box','off','xcolor','none','ycolor','none')
r = 1.1;
for I = 1:13-1
    c = (I-1)/6;
    [num,den] = rat(c);
    tt = c*pi;
    line([0 cos(tt)],[0,sin(tt)],'color',ones(1,3)*0.8)
    
    if num ~= 1
        str = [num2str(num) '\pi'];
    else
        str = '\pi';
    end
    
    if num == 0
        text(r*cos(tt),r*sin(tt),['$0$'],...
            'interpreter','none','horizontalalignment','center','verticalalignment','middle')
    elseif den == 1
        text(r*cos(tt),r*sin(tt),['$' str '$'],...
            'interpreter','none','horizontalalignment','center','verticalalignment','middle')
    else
        text(r*cos(tt),r*sin(tt),['$\frac{' str '}{' num2str(den) '}$'],...
            'interpreter','none','horizontalalignment','center','verticalalignment','middle')
    end
end

N = find(t>1,1);
v(:,1) = f(x(2,:)');
v(:,2) = f(x(N,:)');
c = get(hp,'Color');
hold on
quiver(x([2 N],1),x([2 N],2),v(1,:)',v(2,:)','color',c{1})
quiver(x([2 N],3),x([2 N],4),v(3,:)',v(4,:)','color',c{2})
hold off
