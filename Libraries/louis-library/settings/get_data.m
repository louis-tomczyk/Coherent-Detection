function [data,shape] = get_data(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : 
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-15
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
%   if NARGIN == 1, then get all rows and columns
%   if NARGIN == 2, then get the selected rows, all columns
%   if NARGIN == 3, then get the selected rows and selected columns
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    filename = varargin{1};

    if nargin == 2
        nrow    = varargin{2};
    elseif nargin == 3
        nrow    = varargin{2};
        ncol    = varargin{3};
    end

    data        = readmatrix(filename);
    shape       = size(data);
    
    if nargin == 2
        data    = data(nrow,:);
    elseif nargin == 3
        data    = data(nrow,ncol);
    end
end