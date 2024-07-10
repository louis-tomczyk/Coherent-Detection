function outs = get_pdapr(what,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_PDAPR - Power / Deviation 
%                     to Average Power Ratio
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-03
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get the PAPR or the DAPR in a signal
%
% ----- INPUTS -----
%   WHAT        (string)    Can be 'papr' or 'dapr'
%   VARARGIN    (arrays)    Signals
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    if strcmp(what,'papr') == 1
        PAPRs = zeros(1,nargin-1);
        for k = 1:nargin-1
            mean_tmp    = mean(varargin{k});
            y_tmp       = max(varargin{k});
            PAPRs(k)    = y_tmp/mean_tmp;
        end
        outs = PAPRs;
    elseif strcmp(what,'dapr') == 1
        DAPRs = zeros(1,nargin-1);
        for k = 1:nargin-1
            mean_tmp    = mean(varargin{k});
            std_tmp     = std(varargin{k});
            DAPRs(k)    = std_tmp/mean_tmp;
        end
        outs = DAPRs;
    end