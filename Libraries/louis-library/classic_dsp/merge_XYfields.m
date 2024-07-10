function E = merge_XYfields(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : MERGE_XYFIELDS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-08-05
%   Version         : 1.0
%
% ----- Main idea -----
%   Replace the PBC function of OPTILUX
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    Ex          = varargin{1};
    Ey          = varargin{2};

    Etot        = Ex;
    Etot.field  = Ex.field+Ey.field;

    if Ex.lambda ~= Ey.lambda
        error('Cannot combine fields with different wavelength')
    end
    

    E.lambda    = Ex.lambda;
    Nsamp       = size(Ex.field,1);
    
    if size(Ex.field,2) == 2 && sum(Ex.field(:,2) == zeros(Nsamp,1)) == Nsamp
        Ex.field(:,2) = [];
    end

    if size(Ey.field,2) == 2 &&  sum(Ey.field(:,1) == zeros(Nsamp,1)) == Nsamp
        Ey.field(:,1) = [];
    end


    E.field(1:Nsamp,1)  = Ex.field;
    E.field(1:Nsamp,2)  = Ey.field;

% 
%     if nargin == 3
%         if varargin{3} == "sep"
%             return
%         else
%             E = rmfield(E,'field');
%             I = (real(Ex.field)+real(Ey.field))/2;%sqrt(2);
%             Q = (imag(Ex.field)+imag(Ey.field))/2;%sqrt(2);
%             E.field = I+1i*Q;
%         end
%     end
% 
%     figure
%     plotfield(E,'--p-',struct('power','tot'))
%     plotfield(Etot,'--p-',struct('power','tot'))