classdef mf
    properties
        f;      %f(1) - function handle, f(2) - mf object 
        mnpq;   %f(1)(X) is a p x q matrix and X is a m x n matrix
        cnt;    %cnt - bounds the number of recursions
        d;      %d - depth of the current mf object
    end
    methods
        function out = mf(varargin)
            if nargin > 0
                if isa(varargin{1},'mf')
                    %copy constructor
                    out.f = varargin{1}.f;
                    out.d = varargin{1}.d;
                    out.mnpq = varargin{1}.mnpq;
                    out.cnt = varargin{1}.cnt;
                elseif nargin > 1 && all(size(varargin{1})==[4 1]) && ...
                        (isa(varargin{2},'function_handle') || ~any(varargin{2}))
                    out.mnpq = varargin{1};
                    out.f{1} = varargin{2};
                    out.d = 0;
                    if nargin > 2 && isa(varargin{3},'mf')
                        out.f{2} = varargin{3};
                        out.d = varargin{3}.d+1;
                    end
                else
                    ME = MException('mf:constructor','Invalid function arguments');
                    throw(ME);
                end
            else
                out.mnpq = ones(4,1);
                out.f{1} = @(x) x; %default function is the identity
                out.d = 0;
            end
        end
        function out = subsref(A,B)
            if ~isa(B,'mf')
                if numel(B.subs) == 2
                    ii = B.subs{1};
                    jj = B.subs{2};
                    if A.d > 0
                        M = numel(ii);
                        N = numel(jj);
                        auxii = eye(A.mnpq(3));
                        auxjj = eye(A.mnpq(4));
                        out = mf([A.mnpq(1:2);M;N],...
                                @(y) auxii(ii,:)*A.f{1}(y)*auxjj(:,jj),...
                                D(auxii(ii,:)*A*auxjj(:,jj)));
                    else
                        M = numel(ii);
                        N = numel(jj);
                        auxii = eye(A.mnpq(3));
                        auxjj = eye(A.mnpq(3));
                        out = mf([A.mnpq(1:2);M;N],...
                                @(y) auxii(ii,:)*A.f{1}(y)*auxjj(:,jj));
                    end
                    return
                else
                    B = B.subs{1};
                end
            end
            if isa(B,'mf')
                if B.d > 0 && A.d > 0
                    out = mf([B.mnpq(1:2);A.mnpq(3:4)],...
                        @(y) A.f{1}(B.f{1}(y)),...
                        subsref(D(A),B)*D(B));
                else
                    out = mf([B.mnpq(1:2);A.mnpq(3:4)],...
                        @(y) A.f{1}(B.f{1}(y)));
                end
            else
                out = A.f{1}(B);
            end
        end
        function out = mtimes(A,B)
            if isa(A,'mf') && isa(B,'mf')
                if A.d > 0 && B.d > 0
                    out = mf(...
                            [A.mnpq(1:3);B.mnpq(4)],...
                            @(x) A.f{1}(x)*B.f{1}(x),...
                            kron(B',speye(A.mnpq(3)))*D(A)+kron(speye(B.mnpq(4)),A)*D(B)...
                            );
                else
                    out = mf(...
                            [A.mnpq(1:3);B.mnpq(4)],...
                            @(x) A.f{1}(x)*B.f{1}(x)...
                            );
                end
            elseif isa(A,'mf')
                if A.d > 0
                    out = mf(...
                            [A.mnpq(1:3);size(B,2)],...
                            @(x) A.f{1}(x)*B,...
                            kron(B',eye(A.mnpq(3)))*D(A));
                else
                    out = mf(...
                            [A.mnpq(1:3);size(B,2)],...
                            @(x) A.f{1}(x)*B);
                end
            elseif isa(B,'mf')
                if B.d > 0
                    out = mf([B.mnpq(1:2);size(A,1);B.mnpq(4)],...
                            @(x) A*B.f{1}(x),...
                            kron(eye(B.mnpq(4)),A)*D(B)...
                            );
                else
                    out = mf([B.mnpq(1:2);size(A,1);B.mnpq(4)],...
                            @(x) A*B.f{1}(x)...
                            );
                end
            end
        end
        function out = D(obj)
            if obj.d > 0
                out = mf(obj.f{2});
            else 
                ME = MException('mf:D','Derivative not defined for mf object');
                throw(ME);
            end
        end
        function out = ctranspose(obj)
            if obj.d > 0
                out = mf(...
                        [obj.mnpq(1:2);obj.mnpq(4);obj.mnpq(3)],...
                        @(x) obj.f{1}(x)',...
                        Kmn(obj.mnpq(3),obj.mnpq(4))*D(obj)...
                        );
            else
                out = mf(...
                        [obj.mnpq(1:2);obj.mnpq(4);obj.mnpq(3)],...
                        @(x) obj.f{1}(x)'...
                        );
            end
            out.cnt = obj.cnt;
        end
        function out = kron(A,B)
            if isa(A,'mf') && isa(B,'mf')
                m = B.mnpq(3);
                n = B.mnpq(4);
                p = A.mnpq(3);
                q = A.mnpq(4);
                if A.d > 0 && B.d > 0
                    out = mf([A.mnpq(1:2);A.mnpq(3:4).*B.mnpq(3:4)],...
                            @(x) kron(A.f{1}(x),B.f{1}(x)),...
                            kron(speye(q),kron(sparse(Kmn(n,p)),speye(m)))*...
                            (kron(vec(A),eye(m*n))*D(B)+kron(eye(p*q),vec(B))*D(A))...
                            );
                else
                    out = mf([A.mnpq(1:2);A.mnpq(3:4).*B.mnpq(3:4)],...
                            @(x) kron(A.f{1}(x),B.f{1}(x))...
                            );
                end
            elseif isa(A,'mf')
                m = size(B,1);
                n = size(B,2);
                p = A.mnpq(3);
                q = A.mnpq(4);
                if A.d > 0
                    out = mf([A.mnpq(1:2);A.mnpq(3:4).*size(B)'],...
                        @(x) kron(A.f{1}(x),B),...
                        kron(speye(q),kron(sparse(Kmn(n,p)),speye(m)))*kron(speye(p*q),sparse(vec(B)))*D(A)...
                        );
                else
                    out = mf([A.mnpq(1:2);A.mnpq(3:4).*size(B)'],...
                        @(x) kron(A.f{1}(x),B)...
                        );
                end
            elseif isa(B,'mf')
                m = B.mnpq(3);
                n = B.mnpq(4);
                p = size(A,1);
                q = size(A,2);
                if numel(B.f) > 1
                    out = mf([B.mnpq(1:2);size(A)'.*B.mnpq(3:4)],...
                        @(x) kron(A,B.f{1}(x)),...
                        kron(speye(q),kron(sparse(Kmn(n,p)),speye(m)))*kron(sparse(vec(A)),sparse(eye(m*n)))*D(B)...
                        );
                else
                    out = mf([B.mnpq(1:2);size(A)'.*B.mnpq(3:4)],...
                        @(x) kron(A,B.f{1}(x))...
                        );
                end
            end
        end
        function out = plus(A,B)
            if isa(A,'mf') && isa(B,'mf')
                if A.d > 0 && B.d > 0
                    out = mf(A.mnpq, @(x) A.f{1}(x)+B.f{1}(x),D(A)+D(B));
                else
                    out = mf(A.mnpq, @(x) A.f{1}(x)+B.f{1}(x));
                end
            elseif isa(A,'mf')
                if A.d > 0
                    out = mf(A.mnpq,@(x) A.f{1}(x)+B,D(A));
                else
                    out = mf(A.mnpq,@(x) A.f{1}(x)+B);
                end
            elseif isa(B,'mf')
                if B.d > 0
                    out = mf(B.mnpq,@(x) B.f{1}(x)+A,D(B));
                else
                    out = mf(B.mnpq,@(x) B.f{1}(x)+A);
                end
            end
        end
        function out = minus(A,B)
            if isa(A,'mf') && isa(B,'mf')
                if A.d > 0 && B.d > 0
                    out = mf(A.mnpq,@(x) A.f{1}(x)-B.f{1}(x),D(A)-D(B));
                else
                    out = mf(A.mnpq,@(x) A.f{1}(x)-B.f{1}(x));
                end
            elseif isa(A,'mf')
                if A.d > 0
                    out = mf(A.mnpq,@(x) A.f{1}(x)-B,D(A));
                else
                    out = mf(A.mnpq,@(x) A.f{1}(x)-B);
                end
            elseif isa(B,'mf')
                if B.d > 0
                    out = mf(B.mnpq,@(x) -B.f{1}(x)+A,-D(B));
                else
                    out = mf(B.mnpq,@(x) -B.f{1}(x)+A);
                end
            end
        end
        function out = uminus(A)
            if isa(A,'mf')
                if A.d > 0
                    out = mf(A.mnpq,@(x) -A.f{1}(x),-D(A));
                else
                    out = mf(A.mnpq,@(x) -A.f{1}(x));
                end
            end
        end
        function out = vec(A)
            if isa(A,'mf')
                if A.d > 0
                    out = mf([A.mnpq(1:2);prod(A.mnpq(3:4));1],...
                        @(x) vec(A.f{1}(x)),...
                        D(A)...
                        );
                else
                    out = mf([A.mnpq(1:2);prod(A.mnpq(3:4));1],...
                        @(x) vec(A.f{1}(x))...
                        );
                end
            else
                out = A(:);
            end
        end
        function out = inv(A)
            if isempty(A.cnt)
                A.cnt = A.d;
            end
            if A.d > 0 && A.cnt
                A.cnt = A.cnt-1;
                out = mf(A.mnpq,...
                    @(x) inv(A.f{1}(x)),...
                    -kron(inv(A'),inv(A))*D(A));
            else
                out = mf(A.mnpq,...
                    @(x) inv(A.f{1}(x)));
            end
        end
        function out = get(obj,p)
            switch p
                case 'mnpq'
                    out = obj.mnpq;
                case 'd'
                    out = obj.d;
                case 'f'
                    out = obj.f{1};
                case 'cnt'
                    out = obj.cnt;
            end
        end
        function out = vertcat(varargin)
            if nargin > 2
                out = vertcat(varargin{1},vertcat(varargin{2:end}));
            elseif nargin == 2 && isa(varargin{1},'mf') && isa(varargin{2},'mf')
                A = varargin{1};
                B = varargin{2};
                if ~all(A.mnpq([1 2 4]) == B.mnpq([1 2 4]))
                    ME = MException('mf:vertcat','Dimensions of objects being concatenated is not consistent');
                    throw(ME);
                end
                if A.d > 0 && B.d > 0
                    out = mf([A.mnpq(1:2);A.mnpq(3)+B.mnpq(3);A.mnpq(4)],...
                            @(x) [A.f{1}(x);B.f{1}(x)],...
                            Kmn(A.mnpq(4),A.mnpq(3)+B.mnpq(3))*vertcat(D(A'),D(B')));
                else
                    out = mf([A.mnpq(1:2);A.mnpq(3)+B.mnpq(3);A.mnpq(4)],...
                        @(x) [A.f{1}(x);B.f{1}(x)]);
                end
            else
                ME = MException('mf:vertcat','All input arguments to vertcat must be mf-objects');
                throw(ME);
            end
        end
        function out = horzcat(varargin)
            if nargin > 2
                out = horzcat(varargin{1},horzcat(varargin{2:end}));
            elseif nargin == 2 && isa(varargin{1},'mf') && isa(varargin{2},'mf')
                A = varargin{1};
                B = varargin{2};
                if ~all(A.mnpq([1 2 3]) == B.mnpq([1 2 3]))
                    ME = MException('mf:horzcat','Dimensions of objects being concatenated is not consistent');
                    throw(ME);
                end
                if A.d > 0 && B.d > 0
                    out = mf([A.mnpq(1:2);A.mnpq(3);A.mnpq(4)+B.mnpq(4)],...
                            @(x) [A.f{1}(x);B.f{1}(x)],vertcat(D(A),D(B)));
                else
                    out = mf([A.mnpq(1:2);A.mnpq(3)+B.mnpq(3);A.mnpq(4)],...
                        @(x) [A.f{1}(x);B.f{1}(x)]);
                end
            else
                ME = MException('mf:horzcat','All input arguments to horzcat must be mf-objects');
                throw(ME);
            end
        end
        function out = reshape(A,sz)
            if size(sz,1) ~= 1 && all(~rem(sz,1))
                ME = MException('mf:reshape','sz must be a row vector with integer elements.');
                throw(ME);
            elseif prod(sz) == prod(A.mnpq(3:4))
                if A.d > 0
                    out = mf([A.mnpq(1:2);sz'], @(x) reshape(A.f{1}(x),sz),D(A));
                else
                    out = mf([A.mnpq(1:2);sz'], @(x) reshape(A.f{1}(x),sz));
                end
            else
                ME = MException('mf:reshape','The output must have the same number of elements as the input.');
                throw(ME);
            end
        end
    end
end