function En = fn(Enn,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FN - FIELD NORMALISATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2022-09-09: creation
%   Date            : 2024-01-19: option to remove 0 offset
%   Version         : 1.1
%
% ----- Main idea -----
%   Normalise the fields by wether their own \sqrt{mean power} or by a
%   given \sqrt{mean power} (OPTIONNAL).
%
% ----- INPUTS -----
%   ENN:    (structure) containing the Fields to be normlised (Not Normalised)
%               - LAMBDA    [nm]        wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    
    % checking number of input arguments
    assert(nargin <4,"Wrong number of input arguments." + ...
        " Should be (ENN) or (ENN,PNORM) and optionnaly 'remove' option")

    % input type flexibility
    if isstruct(Enn) == 0
        tmp.lambda = 1550;
        tmp.field  = Enn;
        clear Enn
        Enn = tmp;
    end
    
    if nargin == 2 && check_class_charstring(varargin{1}) == false
        Pnorm = varargin{1};
    end
     
    En.lambda = Enn.lambda;
    
    % n_polar = 2;
    if size(Enn.field,2) == 2*size(Enn.lambda,1)
    
        % split the polarisations
        [Ennx,Enny]= sep_XYfields(Enn);
        Enx = Ennx;
        Eny = Enny;
            
        % normalise each polarisation
        Xnn = Ennx.field;
        Ynn = Enny.field;

        if nargin == 1
            Xn = Xnn/sqrt(mean(empower(Xnn)));
            Yn = Ynn/sqrt(mean(empower(Ynn)));
        elseif class(varargin{1}) == "double"
            Xn = Xnn/sqrt(Pnorm);
            Yn = Ynn/sqrt(Pnorm);
        elseif check_class_charstring(varargin{1}) == true
            Xn = Xnn;
            Yn = Ynn;
        end

        Enx.field = Xn;
        Eny.field = Yn;

        % recombine both polarisations
        En    = merge_XYfields(Enx,Eny);
        
    % n_polar = 1;
    else
        if nargin == 1
            En.field = Enn.field/sqrt(mean(empower(Enn)));
        else
            En.field = Enn.field/sqrt(Pnorm);
        end
    end

    if nargin == 2 && check_class_charstring(varargin{1}) == true
        if strcmp(varargin{1},'remove')== 1

            tmpx = En.field(:,1);
            tmpy = En.field(:,2);

            tmpx = tmpx(tmpx~=0).';
            tmpy = tmpy(tmpy~=0).';

            clear En
            En.lambda       = Enn.lambda;
            En.field(:,1)   = tmpx;
            En.field(:,2)   = tmpy;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   flag = check_class_charstring(in)
% ---------------------------------------------

function flag = check_class_charstring(in)

    if class(in) == "char" || class(in) == "string"
        flag = true;
    else
        flag = false;
    end
