function out = get_value_from_filename(folderPath,quantity,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_VALUE_FROM_FILENAME
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-06
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get a value contained in the name of a file
%
% ----- INPUTS -----
%   FOLDERPATH  (string)    Folder that contain the files
%   QUANTITY    (string)    Keyword before/after which is located the figure
%   VARARGIN    (scalars)[optiona]
%               {1}          if given should be N_CHAR_BEFORE
%               {2}          if given should be N_CHAR_AFTER
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
    
    PathInit = '~/Documents/6_Telecom_Paris/3_Codes/louis/Optilux/';
    cd(PathInit)
    cd(folderPath)
    my_addpath_to_libraries(folderPath,'Optilux',"Libraries/louis library/settings/")
    
    nfiles          = length(dir(pwd))-2;
    folder_struct   = dir(pwd);
    
    
    out = zeros(nfiles,1);
    
    if nargin >= 3
        n_char_b4 = varargin{1};
        if nargin == 4
            n_char_after = varargin{2};
        else
        n_char_after    = 6;
        end
    else
        n_char_b4       = 6;
        n_char_after    = 6;
    end
    
    for k=1:nfiles
    
        filename = folder_struct(k+2).name;
        out(k) = get_number_from_string(filename,quantity,n_char_b4,n_char_after);
    
    end