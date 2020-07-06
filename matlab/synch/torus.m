N1 = 25;
N2 = 25;
N3 = 40;
tt1 = linspace(0,2*pi,N1);
tt2 = linspace(0,2*pi,N2);
tt3 = linspace(0,2*pi,N3);

X = NaN(N2,N1,3);
Y = NaN(N2,N1,3);
Z = NaN(N2,N1,3);
U = NaN(N2,N1,2);
V = NaN(N2,N1,2);
W = NaN(N2,N1,2);
B = NaN(N2,1,2);

for I = 1:N1
    for J = 1:N2
        [x,y,z] = torus_chart(tt1(I),tt2(J));
        x1 = [cos(tt1(I));sin(tt1(I))];
        x2 = [cos(tt1(I)+tt2(J));sin(tt1(I)+tt2(J))];
        for K = 1:2
            q  = 2*K-3;
            xi = [x1;x2;zeros(4,1);q];
            if C(xi)
                aux= law(xi);
                vx = aux(1)+aux(3)*cos(tt1(I));
                vy = aux(2)+aux(3)*sin(tt1(I));
                vz = aux(4);
                U(J,I,K) = vx/norm([vx;vy;vz]);
                V(J,I,K) = vy/norm([vx;vy;vz]);
                W(J,I,K) = vz/norm([vx;vy;vz]);
                X(J,I,K) = x;
                Y(J,I,K) = y;
                Z(J,I,K) = z;
            end
        end
        if ~isnan(X(J,I,1)) && ~isnan(X(J,I,2))
            X(J,I,3) = x;
            Y(J,I,3) = y;
            Z(J,I,3) = z;
        end
    end
end

colors = get(gca,'colororder');
for I=1:3
    hs(I) = surf(X(:,:,I),Y(:,:,I),Z(:,:,I),'facecolor','w','edgecolor',colors(I,:),'linewidth',I);
    hold on
    if I < 3
        quiver3(X(:,:,I),Y(:,:,I),Z(:,:,I),U(:,:,I),V(:,:,I),W(:,:,I),'color',colors(I,:),'linewidth',2)
    end
end
hs(4) = plot3(cos(tt3)+cos(tt3).^2,sin(tt3)+cos(tt3).*sin(tt3),sin(tt3),'k','linewidth',2);
hs(5) = plot3(cos(tt3)+cos(tt3+pi).*cos(tt3),sin(tt3)+cos(tt3+pi).*sin(tt3),sin(tt3+pi),'k--','linewidth',2);
hold off;
axis equal
set(gca,'visible','off','position',[0 0 1 1])
set(gcf,'color','w')
[hl,p1,p2,p3] = legend(hs,{'$\qs=-1$','$\qs=1$','$\qs=\pm 1$','$x=y$','$x=-y$'},...
    'position',[0.6054    0.8214    0.3679    0.1619]);