function out = get_number_from_string(stringIn,what,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_NUMBER_FROM_STRING
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-05
%   Version         : 1.2
%
% ----- MAIN IDEA -----
%   Find in a string the number before a substring.
%
% ----- INPUTS -----
%   stringIn    (string)    string in which we are looking for the number
%   what        (string)    substring located just after the number
%   varargin    (scalar)[optional]
%                           if given, it is the range of characters we extend the search
%
% ----- BIBLIOGRAPHY -----
%   Functions           :
%   Author              : Cedric WANNAZ
%   Author contact      : Contact sheet on Mathworks
%   Date                : 2014-06-20
%   Title of program    : NA
%   Code version        : 1
%   Type                : Matlab Answers
%   Web Address         : https://fr.mathworks.com/matlabcentral/answers/136724-how-do-i-call-the-first-number-in-a-string
% ----------------------------------------------

    iwhat       = strfind(stringIn,what);
    if nargin == 2
        strTmp  = stringIn(iwhat:iwhat+10);
    else
        if nargin > 2
            if iwhat-varargin{1}<1
                istart = 1;
            else
                istart = iwhat-varargin{1};
            end
            if nargin == 4
                if iwhat+varargin{2}>length(stringIn)
                    iend = length(stringIn);
                else
                    iend = iwhat+varargin{2};
                end
            end
            strTmp  = stringIn(istart:iend);
        end
    end

    indexes = regexp(strTmp,'[0123456789.]');
    if sum(indexes-length(strTmp)==-3)==1
        indexes = indexes(1:end-1);
    end

    out = str2double(strTmp(indexes));









