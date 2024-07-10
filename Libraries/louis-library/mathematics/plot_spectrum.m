function plot_spectrum(E,varargin)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : plot_spectrum
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2024-01-16: creation
%   Date            : 2024-01-23
%   Version         : 1.2
%
% ----- MAIN IDEA -----
%   Plot the amplitude spectrum of a field
%
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    if isstruct(E) == 0
        tmp         = E;
        clear E
        E.field     = tmp;
        E.lambda    = 1550;
    end

    x       = E.field;
    xft     = FFT(x);
    xft2    = empower(xft);
    xft2dB  = 10*log10(xft2);
    my_title= inputname(1);

    hold on
    if nargin == 1
        semilogy(xft2dB)
        ylabel("$|FT|^2 [dB/sample]$")
        title(my_title)
        set(gca,"fontsize",15,"FontWeight","bold")
    
    elseif nargin == 2
        if isstruct(varargin{1}) == 0
            if length(varargin{1}) == 2
                ylim(varargin{1})
                semilogy(xft2dB)
            else
                freq = varargin{1};
                plot(freq,xft2dB)
            end
        else
            freq = varargin{1}.freq;
            if freq(end)>1e9
                freq = freq*1e-9;
            end
            semilogy(freq,xft2dB)
        end

        ylim([min(min(xft2dB))-5,max(max(xft2dB))+5])
        xlabel("frequency [GHz]")
        ylabel("$|FT|^2 [dB/GHz]$")

        title(my_title)
        set(gca,"fontsize",15,"FontWeight","bold")

    elseif nargin == 3
        if length(varargin{1}) == 2
            ylim(varargin{1})
            freq = varargin{2}.freq;
        else
            ylim(varargin{2})
            freq = varargin{1}.freq;
        end

        if freq(end)>1e9
            freq = freq*1e-9;

        end

        semilogy(freq,xft2dB)
        xlabel("frequency [GHz]")
        ylabel("$|FT|^2 [dB/GHz]$")
        title(my_title)

        set(gca,"fontsize",15,"FontWeight","bold")
    end

end


%% TRASH CODE

% semilogy(freq,smoothdata(xft2dB,'gaussian',150))