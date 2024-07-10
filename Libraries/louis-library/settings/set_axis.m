function Axis = set_axis(tx,varargin)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_AXIS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   ArXivs          : 2023-01-03: to get only integer Nsps
%   Date            : 2024-01-26: managing SI units
%   Version         : 1.5
%
% ----- Main idea -----
%   Set the Axis structure from basic axis parameters.
%
% ----- INPUTS -----
%   AXIS:   (structure) should be organised as:
%               - NSYMB     (scalar)[]      Number of symbols to send
%               - SYMBRATE  (scalar)[GBAUD] Symbol rate
%
% ----- OUTPUTS -----
%   AXIS:   (structure) to which will be added
%               - SYMBRATE  (scalar)[GBAUD]    Symbol rate
%               - BITRATE   (scalar)[GBIT/S]   Bit rate
%               - FMAX      (scalar)[GHZ]      Maximum frequency
%               - FS        (scalar)[GHZ]      Sampling rate
%               - CS        (scalar)[]         Shannon Coefficient
%               - NSYMB     (scalar)[]         Number of symbols
%               - NBITS     (scalar)[]         Number of bits
%               - NBPS      (scalar)[]         Number of Bits Per Symbol
%               - NSPB      (scalar)[]         Number of Samples Per Bit
%               - NSPS      (scalar)[]         Number of Samples Per Symbol
%               - NSAMP     (scalar)[]         Total number of samples
%               - TSYMB     (scalar)[ns]       Symbol duration
%               - TBIT      (scalar)[ns]       Bit duration
%               - DF        (scalar)[GHZ]      Frequency resolution
%               - DT        (scalar)[ns]       Time resolution
%               - FREQ      (array)[GHZ]       Frequency axis
%               - TIME      (array)[ns]        Time Axis
%
% ----- BIBLIOGRAPHY -----
%   Articles/Books
%   Authors             : louis tomczyk
%   Title               : Transformée de Fourier
%   Publisher           :
%   Volume - N°/edition :
%   Date                : 
%   DOI/ISBN            :
% ---------------------------------------------

    if nargin == 2
        Axis = varargin{1};
    else
        assert(nargin == 1,"nargin == 1 : TX, nargin == 2 : TX,AXIS")
        Axis = struct();
    end
    
    % ==== ROLL OFF
    if isfield(Axis,'roll_off') == 0   
        rolloff = 0.2; % roll off
    else
        rolloff = Axis.rolloff;
    end

    % === TIMES
    if isfield(Axis,'symbrate') == 0
        Axis.symbrate = 40;                     % [GBd]
    end

    Axis.Tsymb  = 1/Axis.symbrate;              % [ns]
    if isfield(Axis,'Nsymb') == 0
        Axis.Nsymb = 2^14;
    end
    Axis.tmax   = Axis.Nsymb*Axis.Tsymb;        % [ns]

    % === FREQUENCIES
    % if even number of channels
    if mod(tx.Nch,2) == 0 
        Axis.fmax   = (tx.Nch-1)/2*tx.deltaf;
    % if odd number of chanels
    else
        fmax_mono   = (1+rolloff)/2*Axis.symbrate;
        Axis.fmax   = 2*floor(tx.Nch/2)*tx.deltaf+fmax_mono;
    end


    % === SHANNON COEFFICIENT
    % ---- Taking into account only Kerr effect
    if isfield(Axis,'cs') == 1
        assert(Axis.cs>=2)
    else        
        Axis.cs  = 6;
    end

    % we want the number of samples per symbol to be even
    % Nsps= cs*fmax/Rs
    Nsps_tmp = Axis.cs*Axis.fmax/Axis.symbrate;

    % if Nsps_tmp not an integer
    if mod(Nsps_tmp,2) ~= 0
        % we look if it ceil value is odd
        ntmp = ceil(Nsps_tmp);
        if mod(ntmp,2) ~= 0
            % if this ceil value is odd, we add 1
            ntmp = ntmp+1;
        end
        % we get our updated fmax
        Axis.fmax= ntmp*Axis.symbrate/Axis.cs;
    end

    % ---- Taking into account Rayleigh, Brillouin, Raman, FWM etc.
    % to come

    Axis.fs     = Axis.cs*Axis.fmax; % [GHz]
    Axis.dt     = 1/Axis.fs;                    % [ns]

    % --- number of Bits Per Symbol
    Axis.Nbps   = tx.Nbps;

    %--- number of bits
    Axis.Nbits  = Axis.Nbps*Axis.Nsymb;

    %--- number of Samples Per Bit
    Axis.Nspb   = Axis.Tsymb*Axis.fs/Axis.Nbps;

    % --- number of Samples Per Symbol
    Axis.Nsps   = Axis.Nbps*Axis.Nspb;

    % --- total number of samples    
    Axis.Nsamp  = Axis.Nsymb*Axis.Nsps;
    if mod(Axis.Nsamp,2)~=0
        if mod(Axis.Nsamp,2)==1
            Axis.Nsamp = Axis.Nsamp+1;
        else
            Axis.Nsamp = ceil(Axis.Nsamp);
        end
    end

    %###% checking that calculations corresponds for the number of samples
    assert(Axis.Nsamp == ceil(Axis.Nsymb*Axis.Nsps))

    % --- frequency resolution
    Axis.df     = Axis.fs/Axis.Nsamp;

    Axis.bitrate= Axis.symbrate*Axis.Nbps;                      % [Gbit/s]
    Axis.Tbit   = Axis.Tsymb/Axis.Nbps;                         % [ns]

    %###% checking that calculations corresponds for the bit time
    assert(Axis.Tbit-Axis.Nspb*Axis.dt<eps)

    %###% time and frequency axis
    Axis.time    = 0:Axis.dt:(Axis.Nsamp-1)*Axis.dt;                % [GHz]
    Axis.freq    = (-Axis.fs/2:Axis.df:Axis.fs/2-Axis.df)+Axis.df/2;% [ns]

    Axis.tmax   = Axis.time(end);
    Axis        = sort_struct_alphabet(Axis);

    if isfield(Axis,"type") == 0
        Axis.type   = "optilux";
    end

    if strcmp(Axis.type,"SI") == 1
        Axis = convert_scales_axis_in(Axis);
    end








function ax = convert_scales_axis_in(Axis,varargin)

    ax      = Axis;
    norm    = 1e9;

    if nargin == 1
        if strcmp(Axis.type,"SI") == 1
            what = "SI";
        else
            what = "optilux";
        end
    else
        what = varargin{1};
    end

    if strcmp(what,"SI") == 1

        ax.symbrate = Axis.symbrate*norm; 
        ax.bitrate  = Axis.bitrate*norm;

        ax.df       = Axis.df*norm;
        ax.fmax     = Axis.fmax*norm;
        ax.fs       = Axis.fs*norm;
        ax.freq     = Axis.freq*norm;

        ax.dt       = Axis.dt/norm;
        ax.tmax     = Axis.tmax/norm;
        ax.Tsymb    = Axis.Tsymb/norm;
        ax.Tbit     = Axis.Tbit/norm;
        ax.time     = Axis.time/norm;

        ax          = sort_struct_alphabet(ax);
        ax.type     = "SI";
        ax.data     = "simu";

    elseif strcmp(what,"optilux") == 1

        ax.symbrate = Axis.symbrate*norm; 
        ax.bitrate  = Axis.bitrate*norm;

        ax.df       = Axis.df*norm;
        ax.fmax     = Axis.fmax*norm;
        ax.fs       = Axis.fs*norm;
        ax.freq     = Axis.freq*norm;

        ax.dt       = Axis.dt/norm;
        ax.tmax     = Axis.tmax/norm;
        ax.Tsymb    = Axis.Tsymb/norm;
        ax.Tbit     = Axis.Tbit/norm;
        ax.time     = Axis.time/norm;

        ax          = sort_struct_alphabet(ax);
        ax.type     = "optilux";
        ax.data     = "simu";
    end




%% TRASH
%     Axis.freq   = linspace(-Axis.fs/2,Axis.fs/2,Axis.Nsamp);    % [GHz]
%     Axis.time   = linspace(0,Axis.tmax,Axis.Nsamp);             % [ns]