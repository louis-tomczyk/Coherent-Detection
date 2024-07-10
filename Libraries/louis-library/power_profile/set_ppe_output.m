function varargout = set_ppe_output(repet_Ntries,method_BP,PP)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_PPE_OUTPUT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-12-01
%   Version         : 1.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
    if repet_Ntries > 1
        % --- back propagation only in the link
        if strcmp(method_BP,'no_pd') == 1
            Mean_PP         = PP{1,2}{1};
            Mean_PP_av      = PP{1,3}{1};

            varargout{1}    = Mean_PP; 
            varargout{2}    = Mean_PP_av;

        else % --- back propagation in link + predisp
            Mean_PP         = PP{1,2}{1};
            Mean_PP_trunc   = PP{1,3}{1};
            Mean_PP_trunc_av= PP{1,4}{1};

            varargout{1}    = Mean_PP; 
            varargout{2}    = Mean_PP_trunc;
            varargout{3}    = Mean_PP_trunc_av;
        end
    else % Ntries = 1
        if strcmp(method_BP,'no_pd') == 1
            PP_av           = PP{1,2}{1};

            varargout{1}    = PP_av;
        else % --- back propagation in link + predisp
            PP_trunc        = PP{1,2}{1};
            PP_trunc_av     = PP{1,3}{1};

            varargout{1}    = PP_trunc;
            varargout{2}    = PP_trunc_av;
        end
    end

end