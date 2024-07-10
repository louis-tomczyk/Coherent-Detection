function is_arg_missing(arg_missing,all_args)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : IS_ARG_MISSING
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-08-10
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Check if argument is missing
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    assert(~isempty(find(all_args == arg_missing, 1)),...
        sprintf('\n     --- %s --- is missing',upper(arg_missing)))

end