function dec = get_decimals(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_DECIMALS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2023-01-31
%   Version         : 2.0
%
% ----- MAIN IDEA -----
% Return the decimal part of a given input
% VARARGIN{1}   the number from which we want the decimals
% VARARGIN{2}{optinal}
%               the number of decimals wanted
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    input   = varargin{1};

    % conversion to string
    if nargin == 1
        str = num2str(input);
    else % choose the number of decimals
        str = num2str(input,varargin{2});
    end

    % locating the decimals
    index   = strfind(str,'.');
    dec.n   = length(char(str(index+1:end)));

    % getting the decimals
    dec.val = str2double(str(index+1:end));

    % if no decimals
    if isnan(dec.val)
        dec.dec = [];
    end
end