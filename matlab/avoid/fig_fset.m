%{ 
function fig_fset(q,Nr,Ntt,rlim,color)
    global delta
    r = linspace(rlim(1),rlim(2),Nr);
    tt = linspace(0,2*pi,Ntt);
    X = zeros(Nr,Ntt);
    Y = zeros(Nr,Ntt);
    Z = zeros(Nr,Ntt);
    for I = 1:Nr
        for J = 1:Ntt
            x = [r(I);cos(tt(J));sin(tt(J))];
            z = fx(x);
            X(I,J) = z(1);
            Y(I,J) = z(2);
            [Vx,minVx] = V(q,x);
            Z(I,J) = Vx-minVx;
        end
    end
    contour(X,Y,Z,[delta delta],'edgecolor',color,'facecolor','none')
    %quiver(X,Y,U,V)
end
%}
function out = fig_fset(q,Nr,Ntt,rlim,color,linestyle,linewidth)
    %grid on
    r = linspace(rlim(1),rlim(2),Nr);
    tt = linspace(0,2*pi,Ntt);
    X = NaN(Nr,Ntt);
    Y = NaN(Nr,Ntt);
    U = NaN(Nr,Ntt);
    V = NaN(Nr,Ntt);
    for I = 1:Nr
        for J = 1:Ntt
            x = [r(I);cos(tt(J));sin(tt(J))];
            xi= [x;q;zeros(2,1)];
            if all(q~=[-1 1]) && C([x;1;zeros(2,1)]) && C([x;-1;zeros(2,1)])
                z = fx(x);
                X(I,J) = z(1);
                Y(I,J) = z(2);
                aux = [zeros(2,4) eye(2)]*f(xi);
                U(I,J) = aux(1)*norm(aux)^(-0.75);
                V(I,J) = aux(2)*norm(aux)^(-0.75);
            elseif any(q==[-1 1]) && C(xi)
                z = fx(x);
                X(I,J) = z(1);
                Y(I,J) = z(2);
                aux = [zeros(2,4) eye(2)]*f(xi);
                U(I,J) = aux(1)*norm(aux)^(-0.75);
                V(I,J) = aux(2)*norm(aux)^(-0.75);
            end
            %hc = draw_circle('center',[1;0],'radius',exp(r(I))+ee);
            %set(hc,'color',0.9*ones(1,3))
            %hl = line(
        end
    end
    out = mesh(X,Y,X.*0,'edgecolor',color,'facecolor','none',...
        'linestyle',linestyle,'linewidth',linewidth)
    %quiver(X,Y,U,V)
end
%