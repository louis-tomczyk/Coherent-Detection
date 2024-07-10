function set_folder_path(mm,iter,what_charac)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_FOLDER_PATH
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-01
%   Version         : 1.0
%
% ----- Main idea -----
%   Set the paths where to save the data for the Power Profile Estimations
%   
% ----- INPUTS -----
%   MM          (scalar)    Parameter value number to change for the PPE
%   ITER        (scalar)    Iteration for the different realisations for
%                           the PPE
%   WHAT_CHARAC (string)    Can be "EPSILON" or 
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if mm == 1
        if iter == 1
            cd ../../ % we are here: Optilux/Test PPE/
            cd(strcat(pwd,'/simulation_parameters/',what_charac))
        else
            cd(strcat(pwd,'Data/simulation_parameters/',what_charac))
        end
    else
        pwd
        if iter == 1
            cd(strcat(pwd,'/Data/simulation_parameters/',what_charac))
        else
            cd(strcat(pwd,'/simulation_parameters/',what_charac))
        end
    end

end