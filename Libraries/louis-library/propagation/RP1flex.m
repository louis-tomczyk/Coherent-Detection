function Efib = RP1flex(Ein,ft,Axis,RPparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : REGULAR PERTURBATION method - order 1
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@mailo.com
%   Date            : 2022-11-23
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
%   Articles/Books
%   Authors             : Armando VANNUCCI, Paolo SERENA, Alberto BONONI
%   Title               : The RP method: a new tool for the iterative
%                         solution of the nonlinear Schrödinger equation
%   Publisher           : JLT
%   Volume - N°/edition : 20 - 7
%   Date                : July 2002
%   DOI/ISBN            : 10.1109/JLT.2002.800376
%
%   Authors             : Takeo SASAI, Etsushi YAMAZAKI, Yoshiaki KISAKA
%   Title               : Performance limit of fibre lonitudinal power
%                         profile estimation methods
%   Publisher           : 
%   Volume - N°/edition : 
%   Date                : 
%   DOI/ISBN            : 
% ----------------------------------------------

% parameters
L       = ft.length;                % [m] length of the span
dz      = RPparams.dz;              % [m] spatial resolution
% P0      = RPparams.P0;              % [mW] input power
gf      = ft.gf*1e-6;               % [1/(mW.m)] fibre nonlinear parameter
Nsteps  = L/dz;
what    = RPparams.what;
loss    = RPparams.loss;
    
% fields initialisation
Enl     = init_field(Ein);          % nonlinear branches
Efib    = init_field(Ein);          % fibre output

tmp = zeros(size(Enl.field));
% linear branch
Elin    = ChromDisp(Ein,ft,Axis,L,loss,"frequency");

% nonlinear branches
if strcmp(what,'RP0') ~= 1
    if strcmp(what,'RP1 sasai') == 1
        loss = 0; % no loss in linear operator, it appears in gf'
        for k = 1:Nsteps-1
            Ecd_part    = ChromDisp(Ein,ft,Axis,k*dz,loss,"time");

            gfp         = gf*P0*exp(-ft.alphaLin*k*dz/2);
            Ekerr_mod   = ModKerr(Ecd_part,gfp,dz,what);
            
            Enl_tmp     = ChromDisp(Ekerr_mod,ft,Axis,L-k*dz,0,"frequency");
            
            Enl.field   = Enl.field+Enl_tmp.field;
        end
    elseif strcmp(what,'RP1 vannucci') == 1
        gf      = 1;
        loss    = 1;
        parfor k = 1:Nsteps-1
            Ecd_part    = ChromDisp(Ein,ft,Axis,k*dz,loss,"time");
            Ekerr_mod   = ModKerr(Ecd_part,gf,k*dz,what);
            Enl_tmp     = ChromDisp(Ekerr_mod,ft,Axis,L-k*dz,loss,"frequency");
            
%             Enl.field   = Enl.field+Enl_tmp.field;
            tmp         = tmp+Enl_tmp.field;
        end
    end
end

Enl.field = tmp;
if strcmp(what,'RP1 vannucci') == 1
    gf          = ft.gf*1e-6;
    Efib_tmp    = Elin.field-1i*gf*dz*Enl.field;
elseif strcmp(what,'RP1 sasai') == 1
    Efib_tmp    = Elin.field+Enl.field;
end

% get back to time-domain
Efib.field = iFastFT(Efib_tmp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Xft = FastFT(X)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FastFT - Fast Fourier Transform
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-17
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- outPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------


Xft = fftshift(fft(fftshift(X)));

function X = iFastFT(Xft)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : IFastFT - Inverse Fast Fourier Transform
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-17
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- outPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------


X = fftshift(ifft(fftshift(Xft)));

% ----------------------------------------------------------------------- %
function Ecd = ChromDisp(Ein,ft,Axis,z,loss,domain)
% Chromatic Dispersion

    beta2   = ft.beta2;             % [s²/rad²/m]
    om      = 2*pi*Axis.freq*1e9;   % [rad.Hz]

    if loss == 0
        H   = transpose(fastexp(-beta2/2*om.^2*z));
    else
        H   = transpose(fastexp(-beta2/2*om.^2*z).*exp(-ft.alphaLin*z/2));
    end

    Ecd     = Ein;
    X       = Ein.field;
    Xfft    = FastFT(X);
    Xfftcd  = Xfft.*H;

    if strcmp(domain,"time") == 1
        Xcd = iFastFT(Xfftcd);
    elseif strcmp(domain,'frequency') == 1
        Xcd = Xfftcd;
    end
    
Ecd.field = Xcd;

% ----------------------------------------------------------------------- %
function Ekerr = ModKerr(Ecd_part,gfp,dz,what)

% MODified KERR operator

    Ekerr   = Ecd_part;

    if strcmp(what,'RP1 vannucci') == 1
        Ekerr.field = empower(Ecd_part).*Ecd_part.field;
    elseif strcmp(what,'RP1 sasai') == 1
        Ekerr.field = -1i*gfp*dz*empower(Ecd_part).*Ecd_part.field;
    end


% ----------------------------------------------------------------------- %
function E = init_field(model_field)

    E       = model_field;
    E.field = zeros(size(E.field));