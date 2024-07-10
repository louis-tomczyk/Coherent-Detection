function Eds = ds(E,NspsIni,NspsFinal,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : DS - DOWNSAMPLING
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.-tomczyk@telecom-paris.com
%   ArXivs          : 2023-01-10: creation
%   Date            : 2024-01-17: exlpliciting varargin (not 0/1 only, also "replace" available)
%   Version         : 1.2.2
%
% ----- Main idea -----
%   Downsample signals contained in structures as in OPTILUX
%
% ----- INPUTS -----
%   E           (structure) containing the field(s) to downsample
%                   - LAMBDA(scalar)[nm]       wavelength
%                   - FIELD (array)[sqrt(mW)]  normalised electric fields
%   NSPSINI     (scalar)    Numbers of elements before downsampling
%   NSPSFINAL   (scalar)    Numbers of elements after downsampling
%   REPLACEMENT (BOLLEAN){optional}
%                  If == 1, then we replace the erase the not
%                  downsampled field
%
% ----- OUTPUTS -----
%  EDS:         (structure) containing the laser field(s) downsampled
%                   - LAMBDA    (scalar)[nm]       wavelength
%                   - FIELD     (array)[sqrt(mW)]  normalised fields
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    assert(NspsFinal<=NspsIni,['Cannot down sample to less than' ...
        'one sample per symbol. Increase Nsps'])
    
    NN      = ceil(NspsIni/NspsFinal);
    offset  = 0;

    if isfield(E,'lambda') == 0
        E.lambda = 1550;
    end
    
    if size(E.field,2) == 2*size(E.lambda,1)
        
        tmp         = E.field;
        [Ex,Ey]     = sep_XYfields(E);
        Edsx.lambda = E.lambda;
        Edsy.lambda = E.lambda;

        Edsx.field  = Ex.field;
        Edsy.field  = Ey.field;

        Edsx.field  = downsample(Ex.field,NN,offset);
        Edsy.field  = downsample(Ey.field,NN,offset);
        Eds         = merge_XYfields(Edsx,Edsy);
        Eds.field_ds= Eds.field;
        Eds.field   = tmp;
    else

        Eds.lambda  = E.lambda;
        Eds.field   = E.field;
        Eds.field_ds= downsample(E.field,NN);
        
    end

    if nargin == 4 && (strcmp(varargin{1},"replace") ==1 || varargin{1} == 1)
        Eds.field   = Eds.field_ds;
        Eds         = rmfield(Eds,'field_ds');
    end

    Eds.Nsps = NspsFinal;
end