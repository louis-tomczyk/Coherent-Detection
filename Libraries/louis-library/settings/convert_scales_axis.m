function ax = convert_scales_axis(Axis,varargin)
%%
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CONVERT_SCALES_AXIS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   ArXivs          : 
%   Date            : 2024-01-26: creation
%   Version         : 1.0
%
% ----- Main idea -----
%   Change the units of the Axis parameters
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    %% maintenance
    ax      = Axis;
    norm    = 1e9;

    %%% managing output scale
    if nargin == 1
        if strcmp(Axis.type,"SI") == 1
            what = "optilux";
        else
            what = "SI";
        end
    else
        what = varargin{1};
    end

    %%% managing is it is TX or ADC
    if isfield(ax,"Tsymb") == 1
        flag_isADC = false;
    else 
        flag_isADC = true;
    end

    %% conversion
    if strcmp(what,"SI") == 1

        ax.df       = Axis.df*norm;
        ax.freq     = Axis.freq*norm;
        ax.fs       = Axis.fs*norm;
        
        ax.dt       = Axis.dt/norm;
        ax.tmax     = Axis.tmax/norm;
        ax.time     = Axis.time/norm;

        ax.symbrate = Axis.symbrate*norm; 

        if isfield(ax,'fmax') == 1
            ax.fmax     = Axis.fmax*norm;
        end
        
        if flag_isADC == true && isfield(ax,'fs_osc') == 1
            ax.fs_osc = ax.fs_osc*norm;
        end
        
        if flag_isADC == false
            ax.bitrate  = Axis.bitrate*norm;
            ax.Tsymb    = Axis.Tsymb/norm;
            ax.Tbit     = Axis.Tbit/norm;
        end

        ax          = sort_struct_alphabet(ax);
        ax.type     = "SI";

    elseif strcmp(what,"optilux") == 1

        ax.df       = Axis.df/norm;
        ax.freq     = Axis.freq/norm;
        ax.fs       = Axis.fs/norm;

        ax.dt       = Axis.dt*norm;
        ax.tmax     = Axis.tmax*norm;
        ax.time     = Axis.time*norm;

        ax.symbrate = Axis.symbrate/norm;

        if isfield(ax,'fmax') == 1
            ax.fmax     = Axis.fmax/norm;
        end

        if flag_isADC == true && isfield(ax,'fs_osc') == 1
            ax.fs_osc = ax.fs_osc/norm;
        end
            
        if flag_isADC == false
            ax.bitrate  = Axis.bitrate/norm;
            ax.Tsymb    = Axis.Tsymb*norm;
            ax.Tbit     = Axis.Tbit*norm;
        end

        ax          = sort_struct_alphabet(ax);
        ax.type     = "optilux";
    end

