function plot_phase(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PLOT_PHASE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-16
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% VARARGIN{1} --- Field
% VARARGIN{2} --- AXIS
% VARARGIN{3} --- 2 (if 2 pol wanted) [OPTIONAL]
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    E       = varargin{1};
    Axis    = varargin{2};
    Nsamp   = size(E.field,1);

    if Nsamp > Axis.Nsymb
        % getting the number of samples per symbol
        Nsps    = Nsamp/Axis.Nsymb;
        assert(is_integer(Nsps) == 1,"Number of samples per symbol should be an integer")

        % downsample to get 1 sample per symbol
        E       = ds(E,1,Nsps,1);
    end

    % *4/pi if QPSK modulation
    hold on
    plot(angle(E.field(:,1))*4/pi)
    if nargin == 3
        plot(angle(E.field(:,2))*4/pi)
    end
    xlim([0,Axis.Nsymb+1])
    xlabel("Symbol number")
    ylabel("Symbol value")
    set(gca,'fontsize',10,'fontweight','bold')
end