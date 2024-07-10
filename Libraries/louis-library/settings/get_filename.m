function filename = get_filename(directory,keywords)

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
%   Find the name of the file located in the DIRECTORY which contains input
%   KEYWORDS as arguments.
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    cd(directory)
    folder_struct       = dir(pwd);
    nfiles              = length(dir(pwd))-2;
    nkeywords           = length(keywords);

    for k = 1:nfiles
        filename_tmp    = folder_struct(k+2).name;
        count_match     = 0;

        for j = 1:nkeywords
            if contains(filename_tmp,keywords{j}) == 1
                count_match     = count_match+1;
                if count_match == nkeywords
                    filename    = filename_tmp;
                    break;
                end
            end
        end
        if count_match == nkeywords
            filename    = filename_tmp;
            break;
        end
    end

end