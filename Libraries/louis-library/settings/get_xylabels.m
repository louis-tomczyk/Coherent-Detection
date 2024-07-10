function labels = get_xylabels(x,y,varargin)

    N = length(x);
    labels = cell(1,N);

    if nargin >2
        assert(nargin == 4," ======= 2 optional arg for X & Y needed ======= ")
        xres = varargin{1};
        yres = varargin{2};
    end

    if nargin > 2
        for k = 1:N
            labels{k} = strcat('(',num2str(round(x(k),xres)),',',num2str(round(y(k),yres)),')');
        end
    else
        for k = 1:N
            labels{k} = strcat('(',num2str(x(k)),',',num2str(y(k)),')');
        end
    end
    