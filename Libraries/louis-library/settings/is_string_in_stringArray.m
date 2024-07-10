function bool = is_string_in_stringArray(stringArray,myString)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name      IS STRING IN STRINGARRAY
%   Author             louis tomczyk
%   Institution        Telecom Paris
%   Email              louis.tomczyk@telecom-paris.fr
%   Date               2023-02-24
%   Version            1.0
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    assert(nargin == 2,"too much or too few arguments")
    assert(class(stringArray)=="char" || class(stringArray)=="string",...
        "arg 1 should be a string")
    
    bool = isempty(find(ismember(stringArray,myString), 1)) == 0;