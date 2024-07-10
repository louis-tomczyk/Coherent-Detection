function structs_packed = pack_structs(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PACK_STRUCTS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-23
%   Version         : 1.0
%
% ----- Main idea -----
%   Pack into one structure several structures
%
% ----- INPUTS -----
%   VARARGIN : [structures] structures to pack
%
% ----- OUTPUTS -----
%  STRUCTS_PACKED : [structure] containing the input structures
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    % pre-allocate a cell-array which will contain the converted-to-string
    % input variable names
    names = cell(nargin,1);

    % put the names of all the input structures in a cell-array
    for k=1:nargin
        tmp       = inputname(k);
        names(k)  = {tmp(2:end)};
    end

    % pack the strctures
    for k=1:nargin
        structs_packed.(names{k}) = varargin{k};
    end
end