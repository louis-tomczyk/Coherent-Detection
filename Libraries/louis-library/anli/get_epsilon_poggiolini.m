function varargout = get_epsilon_poggiolini(ft,Axis,tx,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_EPSILON_POGGIOLINI
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-02-14
%   Version         : 1.0
%
% ----- Main idea -----
% ----- INPUTS -----
%   FT      (structure)    
%   AXIS    (structure)
%   TX      (structure)
%   VARARGIN(string)        Can be "approx" or nothing
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    Leffa   = 1/2/ft.alphaLin;
    Lspan   = ft.length;
    BW      = (1+tx.rolloff)*Axis.symbrate*1e9;

    if nargin == 4
        if strcmp(varargin{1},'approx') == 1
            varargout{1}= 18/5/(pi^2)/abs(ft.beta2)/Lspan/(BW^2);
            varargout{2}= get_epsilon_poggiolini(ft,Axis,tx);
            varargout{3}= abs(varargout{1}-varargout{2})/varargout{2}*100;
        end
    else
        InAsinh         = pi^2/2*abs(ft.beta2)*Leffa*BW^2;
        Denum           = Lspan*asinh(InAsinh);
        Num             = 6*Leffa;
        InLog           = 1+Num/Denum;

        varargout{1}    = 3/10*log(InLog);
    end

