function DSP = set_dsp(varargin)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_DSP - Digital Signal Processing
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-09-02
%   Version         : 1.3
%
% ----- Main idea -----
%   Set the DSP structure from the laser transmitter parameters and
%   optional dsp parameters
% 
% ----- INPUTS -----
%   DSP:        (strcture) containing DSP parameters - OPTIONAL
%                   - WR        (bool)[]    = 1 if decision step
%                   - CPC_AVG   (scalar)[]  window averaging length
% ----- OUTPUTS -----
%   DSP         (structure) containing the fields:
%               - CDC       (bool)  Chromatic Dispersion Compensation
%               - PMDC      (bool)  Polarisation Mode DC
%               - CPC       (bool)  Carrier Phase C
%               - CPC_AVG   (scalar)Window averaging length
%               - WR        (bool)  Waveform Reconstruction
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    argnames = string(nargin);
    for k=1:nargin
        argnames(k)  = inputname(k);
    end

    DSP.cdc = 1;
    DSP.cpc = 1;
    DSP.mimo    = 1;
    
    if nargin == 0
        DSP.wr      = 0;
        DSP.nlne    = 0;
    end

    if nargin == 1
        is_arg_missing('dsp',argnames);
        dsp = varargin{argnames == 'dsp'};

        if isfield(dsp,'wr') == 0
            DSP.wr = 0;
        else
            DSP.wr = dsp.wr;
        end

        if isfield(dsp,'nlne') == 0
            DSP.nlne = 0;
        else
            DSP.nlne = dsp.nlne;
        end
    end

end