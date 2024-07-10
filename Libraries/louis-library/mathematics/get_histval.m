function [x,y] = get_histval(in,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_HISTVAL
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-08-01
%   Version         : 1.0
%
% ----- Main idea -----
%   Get the curved obtained by the center of each rectangle of the
%   histogram
%
% ----- INPUTS -----
%   INPUT (array)
%   NBINS (scalar){optional}    Length of the X-axis
%
% ----- OUTPUTS -----
%   X   (array)     X values
%   Y   (scalar)    Y values
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    

    assert(nargin >= 1,"missing the vector of values")
    
    if nargin == 1
        H       = histogram(in,'Normalization','count','Visible','off');
    else
        nbins   = varargin{1};
        H       = histogram(in,nbins,'Normalization','count','Visible','off');
    end

    y       = H.Values;
    BinEdges= H.BinEdges;
    x       = (BinEdges(1:end-1)+BinEdges(2:end))/2;
    close all
