function Efilt = filt_opt(Ein,Axis,rx)
%%
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FILT_OPT
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
%%

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain
    % checking power levels
    
    Hopt        = opt_filter(Axis,rx);
    Efiltft     = Hopt.*FFT(Ein,"array");
    Efilt       = iFFT(Efiltft,'struct');
    Efilt.Nsps  = Ein.Nsps;
    
    Diff_pows   = sum(abs(get_power(Ein)-get_power(Efilt)));
    ErrRel_pows = Diff_pows/sum(get_power(Ein));

    assert(ErrRel_pows<5e-2,"power mismatch")





function H_opt = opt_filter(Axis,rx)

    n_polar     = rx.CMA.txpolars;

    if abs(Axis.freq(1))>1e5
        scale = 1e-9;
    else
        scale = 1;
    end

    Hopt    = exp(-(Axis.freq*scale/rx.OEF.opt_bw).^(2*rx.OEF.opt_order));
    H_opt   = repmat(Hopt,[n_polar,1]).';
%-----------------------------------------------------
