function export_structure(varargin)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : EXPORT_STRUCTURE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-09-20
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Export structure content into XML files (XML is imposed by the
%   WRITESTRUCTURE Matlab's function).
%   It is done in a recurssive way as a structure can contain nested
%   sub-structures and the WRITESTRCTURE function does not export those.
%
% ----- INPUTS -----
%   INPUT STRUCTURE (structure) To be exported
%   FILENAME        (string)    Name of the exported file
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    input_structure = varargin{1};
    if nargin == 2
        export_name = varargin{2};
    end
    assert(isstruct(input_structure)==1,"Input is NOT a structure.")

    % 0/ input structure informations
    struct_name = inputname(1);

    % 1/ export all the structures
    if nargin == 1
        writestruct(input_structure, ...
            strcat(string(datetime),sprintf(...
            " --- %s parameters.xml",struct_name)));
    else
        writestruct(input_structure, ...
            strcat(string(datetime),sprintf(...
            " --- %s.xml",export_name)));
    end