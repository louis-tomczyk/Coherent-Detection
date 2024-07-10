function [Ecd,ft] = CD(Ein,Axis,ft)
%%
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CD - CHROMATIC DISPERSION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 
%   Date            : 2024-01-30: creation
%   Version         : 1.0
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    % filter definition in frequency domain
    % filtering in frequency domain
    % back to the time domain

    n_polar = size(Ein.field,2);
    ft.H_CD = CD_filter_in(Axis,ft,n_polar);
    Ecdft   = ft.H_CD.*FFT(Ein,"array");
    Ecd     = iFFT(Ecdft,'struct');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------
% ----- CONTENTS -----
%   H_CD = CD_filter_in(Axis,ft,n_polar)
% ---------------------------------------------

function H_CD = CD_filter_in(Axis,ft,n_polar)
    
    f       = Axis.freq*1e9;
    b2      = ft.beta2;
    L       = ft.length;

    H_CD    = exp(1i*b2/2*(2*pi*f).^2*L);
    H_CD    = repmat(H_CD,[n_polar,1]).';
% ---------------------------------------------