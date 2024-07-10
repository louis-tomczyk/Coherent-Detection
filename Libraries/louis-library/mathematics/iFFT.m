function X = iFFT(Xft,varargin)

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

if isstruct(Xft) == 1
    flag    = 1;
    lambda  = Xft.lambda;
    tmp     = Xft.field;
    clear Xft
    Xft     = tmp;
    clear tmp
else
    lambda  = 1550;
    flag    = 0;
end

% performing the fft

SIZE    = size(Xft);
Nsamp   = max(SIZE);
Npolar  = min(SIZE);

if size(Xft,1) == 2
    Xft = Xft.';
end

Norm    = Nsamp;
X       = fftshift(ifft(fftshift(Xft.*Norm)));

% checking Parseval-Plancherel theorem

power_freq = sum(abs(Xft).^2);
power_time = sum(abs(X).^2)/Nsamp;

assert(sum(round(power_time,3)==round(power_freq,3))==Npolar, ...
    sprintf("Parseval-Plancherel theorem not verified: \n" + ...
    "power time = %.3f [mW]\npower_freq = %.3f [mW]",...
    sum(power_time),sum(power_freq)))


% managing output type

if nargin == 2 && strcmp(varargin{1},"struct") == 1
    tmp         = X;
    clear X
    X.lambda  = lambda;
    X.field   = tmp;
elseif nargin == 2 && strcmp(varargin{1},"struct") == 0
    return
else
    if flag == 1
        tmp         = X;
        clear X
        X.lambda  = lambda;
        X.field   = tmp;
    end
end

