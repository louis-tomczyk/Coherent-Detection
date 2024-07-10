function my_addpath_to_libraries(myPathInit,pathKeyword,pathToAdd)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : MY_ADDPATH__TO_LIBRARIES
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-03-05
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Enable, from any subfolder of the simulator, to find the path to the libraries
%
% ----- INPUTS -----
%   MYPATHINIT  (string)    Where we are currently
%   PATHKEYWORD (string)    Folder name from which we need to count the number of times we need to
%                           mount in the tree
%   PATHTOADD   (string)    Path to add
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    pathTmp     = myPathInit;
    indexPath   = strfind(pathTmp,pathKeyword);
    strTmp      = myPathInit(indexPath:end);
    nFolders    = sum(regexp(strTmp,'[/]')~=0);

    strTmp = string();
    for k=1:nFolders
        strTmp = strTmp+'../';
    end
    pathSettings = strTmp+pathToAdd;
    
    cd(myPathInit)
    addpath(strTmp)
    addpath(pathSettings)