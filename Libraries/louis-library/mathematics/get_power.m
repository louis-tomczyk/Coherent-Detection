function P = get_power(varargin)

%%
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_POWER
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-09-03
%   Version         : 1.5
%
% ----- Main idea -----
%   Calculate the power contained in a field
%
% ----- INPUTS -----
%   VARARGIN{1}:(structure) containing the Field - MANDATORY
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
%   VARARGIN{2}:(string) the power unit - OPTIONAL
%               Possible choices are: mW - W - dBm
%
% ----- OUTPUTS -----
%   P:  (array)[VARARGIN{2}] containing the powers
%               - P(1) the total power, whether it is Polarization
%               Multiplexed or not - SYSTEMATIC
%               - P(2) power in polarization X - OPTIONAL
%               - P(3) power in polarization Y - OPTIONAL
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    
    if isstruct(varargin{1}) == 0
        tmp_f.field = varargin{1};
    else
        tmp_f = varargin{1};
    end

    if size(tmp_f.field,2) == 1
        P = mean(empower(tmp_f));
    else
        P(2) = mean(empower(tmp_f.field(:,1)));
        P(3) = mean(empower(tmp_f.field(:,2)));
        P(1) = P(2)+P(3);
    end

    if nargin == 1
        return
    else
        if isfield(varargin{2},'unit') == 1
            if strcmp(varargin{2}.unit,'W')
                P = P*1e-3;
            elseif strcmp(varargin{2}.unit,'dBm')
                P = 10*log10(P);
            elseif ~strcmp(varargin{2}.unit,'mW')
                tmp         = input("Only [W] or [dBm] are available. [mW] is by default.   >> ",'s');
                options     = varargin{2};
                options.unit= tmp;
                P           = get_power(varargin{1},struct('unit',tmp));
            end
        end

        if isfield(varargin{2},'polar') == 1
            if strcmp(varargin{2}.polar,'x')
                P = P(2);
            elseif strcmp(varargin{2}.polar,'y')
                P = P(3);
            elseif strcmp(varargin{2}.polar,'xy')
                P = [P(2),P(3)];
            elseif strcmp(varargin{2}.polar,'tot')
                P = P(1);
            elseif strcmp(varargin{2}.polar,'all')
                return
            else
                tmp = input("Only [x], [y], [xy] or [tot] are available. [tot] is by default.   >> ",'s');
                P   = get_power(varargin{1},struct('polar',tmp));
            end
        end

        if isfield(varargin{2},'IQ') == 1
            if strcmp(varargin{2}.IQ,'I')
                P = mean(empower(real(tmp_f.field)));

            elseif strcmp(varargin{2}.IQ,'Q')
                P = mean(empower(imag(tmp_f.field)));

            elseif strcmp(varargin{2}.IQ,'XI')
                P = mean(empower(real(tmp_f.field(:,1))));
            
            elseif strcmp(varargin{2}.IQ,'XQ')
                P = mean(empower(imag(tmp_f.field(:,1))));
            
            elseif strcmp(varargin{2}.IQ,'YI')
                P = mean(empower(real(tmp_f.field(:,2))));
            
            elseif strcmp(varargin{2}.IQ,'YQ')
                P = mean(empower(imag(tmp_f.field(:,2))));
            
            elseif strcmp(varargin{2}.IQ,'all_X')
                PXI = mean(empower(real(tmp_f.field(:,1))));
                PXQ = mean(empower(imag(tmp_f.field(:,1))));
                PX  = PXI + PXQ;

                P = [PXI,PXQ,PX];  
            elseif strcmp(varargin{2}.IQ,'all_Y')
                PYI = mean(empower(real(tmp_f.field(:,2))));
                PYQ = mean(empower(imag(tmp_f.field(:,2))));
                PY  = PYI + PYQ;

                P = [PYI,PYQ,PY]; 
            elseif strcmp(varargin{2}.IQ,'all_I')
                PXI = mean(empower(real(tmp_f.field(:,1))));
                PYI = mean(empower(real(tmp_f.field(:,2))));
                PI  = PXI + PYI;

                P = [PXI,PYI,PI]; 
            elseif strcmp(varargin{2}.IQ,'all_Q')
                PXQ = mean(empower(imag(tmp_f.field(:,1))));
                PYQ = mean(empower(imag(tmp_f.field(:,2))));
                PQ  = PXQ + PYQ;

                P = [PXQ,PYQ,PQ]; 
            elseif strcmp(varargin{2}.IQ,'all')
                PXI = mean(empower(real(tmp_f.field(:,1))));
                PXQ = mean(empower(imag(tmp_f.field(:,1))));
                PYI = mean(empower(real(tmp_f.field(:,2))));
                PYQ = mean(empower(imag(tmp_f.field(:,2))));

                PX  = PXI + PXQ;
                PY  = PYI + PYQ;
                PI  = PXI + PYI;
                PQ  = PXQ + PYQ;

                Pt  = PX + PY;
                P = [PXI,PXQ,PYI,PYQ,PX,PY,PI,PQ,Pt];
            else
                tmp = input("Only [I], [Q], [XI], [XQ], [YI], [YQ], [all_I], [all_Q], [all_X], [all_Y], [all] are available.   >> ",'s');
                P   = get_power(varargin{1},struct('IQ',tmp));
            end
        end

    end
    
end