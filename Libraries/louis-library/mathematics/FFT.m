function Xft = FFT(X,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_DISTANCE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2022-09-17: creation
%                   : 2024-01-16: adding normalisation 
%                       and Parseval-Plancherel check
%   Date            : 2024-01-22: adding flexibility of input
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

if isstruct(X) == 1
    flag    = 1;
    lambda  = X.lambda;
    tmp     = X.field;
    clear X
    X       = tmp;
    clear tmp
else
    flag = 0;
end

% performing the fft

SIZE    = size(X);
Nsamp   = max(SIZE);
Npolar  = min(SIZE);

if size(X,1) == 2
    X = X.';
end

Norm    = 1/Nsamp;
Xft     = fftshift(fft(fftshift(X.*Norm)));

% checking Parseval-Plancherel theorem

power_time = sum(abs(X).^2)/Nsamp;
power_freq = sum(abs(Xft).^2);

assert(sum(round(power_time,3)==round(power_freq,3))==Npolar, ...
    sprintf("Parseval-Plancherel theorem not verified: \n" + ...
    "power time = %.3f [mW]\npower_freq = %.3f [mW]",...
    sum(power_time),sum(power_freq)))


% managing output type

if nargin == 2 && strcmp(varargin{1},"struct") == 1
    tmp         = Xft;
    clear Xft
    Xft.lambda  = lambda;
    Xft.field   = tmp;
elseif nargin == 2 && strcmp(varargin{1},"struct") == 0
    return
else
    if flag == 1
        tmp         = Xft;
        clear Xft
        Xft.lambda  = lambda;
        Xft.field   = tmp;
    end
end
