function varargout = get_tofolder(getting_data,stdized,what_charac,Alphaa)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_TOFOLDER
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-15
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Whether go to the folder where are stored the data we want to study,
%   typically the MEANS_AND_STDS folder.
%   Whether go to the folder where the raw data are to make operations on
%   it, typically getting the basic statistics on those.
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    alphaa = strcat("alpha_",num2str(Alphaa));

    if getting_data == 1
        if stdized == 0
            if strcmp(what_charac,"epsilon") == 1
    
                cd(strcat("Data/",alphaa,"/txt/not stdized/",what_charac,'/'))
                varargout{2}  = length(dir(pwd))-2; % -2 because of  "./" and "../" files
                assert(varargout{2} ~= 0,"No files in this directory")
                varargout{3} = zeros(varargout{2},1);
    
            elseif strcmp(what_charac,"pd") == 1
    
                cd(strcat("Data/",alphaa,"/txt/not stdized/",what_charac,'/'))
                varargout{2}    = length(dir(pwd))-2;
                assert(varargout{2} ~= 0,"No files in this directory")
                varargout{3}    = zeros(varargout{2},1);
    
            end
        else
            if strcmp(what_charac,"epsilon") == 1
    
                cd(strcat("Data/",alphaa,"/txt/stdized/",what_charac,'/'))
                varargout{2}  = length(dir(pwd))-2;
                assert(varargout{2} ~= 0,"No files in this directory")
                varargout{3} = zeros(varargout{2},1);
    
            elseif strcmp(what_charac,"pd") == 1
    
                cd(strcat("Data/",alphaa,"/txt/stdized/",what_charac,'/'))
                varargout{2}  = length(dir(pwd))-2;
                assert(varargout{2} ~= 0,"No files in this directory")
                varargout{3} = zeros(varargout{2},1);
    
            end
        end
    else % if getting data == 0
        if stdized == 1
            if strcmp(what_charac,"epsilon") == 1
                cd(strcat(pwd,"/Data/",alphaa,"/txt/stdized/means_and_stds/",what_charac))
            elseif strcmp(what_charac,"pd") == 1
                cd(strcat(pwd,"/Data/",alphaa,"/txt/stdized/means_and_stds/",what_charac))
            end
        else
            if strcmp(what_charac,"epsilon") == 1
                cd(strcat(pwd,"/Data/",alphaa,"/txt/not stdized/means_and_stds/",what_charac))
            elseif strcmp(what_charac,"pd") == 1
                cd(strcat(pwd,"/Data/",alphaa,"/txt/not stdized/means_and_stds/",what_charac))
            end
        end
    end

    varargout{1} = dir(pwd);

end