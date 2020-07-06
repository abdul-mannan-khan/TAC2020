function hdl = create_axis(layout,width,varargin)
    top = 0;
    bot = 0.1;
    lft = 0.1;
    rgt = 0;
    inx = 0.025;
    iny = 0.025;
    for I = 1:2:numel(varargin)
        if strcmpi(varargin{I},'TopMargin')
            top = varargin{I+1};
        end
        if strcmpi(varargin{I},'BottomMargin')
            bot = varargin{I+1};
        end
        if strcmpi(varargin{I},'LeftMargin')
            lft = varargin{I+1};
        end
        if strcmpi(varargin{I},'RightMargin')
            rgt = varargin{I+1};
        end
        if strcmpi(varargin{I},'InnerXMargin')
            inx = varargin{I+1};
        end
        if strcmpi(varargin{I},'InnerYMargin')
            iny = varargin{I+1};
        end
    end
    [m,n] = size(layout);
    M = max(max(layout));
    figure('units','centimeters','position',[0 0 width m*width/n])
    for I = 1:M
        [i,j] = find(layout == I);
        i0 = min(i);
        j0 = min(j);
        i1 = max(i);
        j1 = max(j);
        w0 = (1-lft-rgt)*(j1-j0+1)/n;
        h0 = (1-top-bot)*(i1-i0+1)/m;
        x = (1-lft-rgt)*(j0-1)/n+lft+inx;
        w = w0-2*inx;
        y = (1-bot-top)*(m-i1)/m+bot+iny;
        h = h0-2*iny;
        hdl(I) = axes('position',[x y w h]);
    end