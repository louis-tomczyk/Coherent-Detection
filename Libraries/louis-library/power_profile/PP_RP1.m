function PP = PP_RP1(Ein,Eoptfilt,ft,amp,RPparams)

% function PP = PP_RP1(Ein,Eoptfilt,ft,amp,RPparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : RP1STEP
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-27
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Kerr step sandwiched by two dispersion steps
%
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
%   Functions   : FIBER
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% -----------------------
%   Articles
%   Author              :
%   Title               :
%   Jounal              :
%   Volume - N°         :
%   Date                :
%   DOI                 :
% ----------------------------------------------

Nsamp   = size(Ein.field,1);
Nsteps  = amp.Nspan*ft.length/RPparams.dz;
PP      = zeros(Nsamp,Nsteps);

for k = 1:Nsteps
    disp("==========================================")

    RPparams.dz_part        = k*RPparams.dz;
    RPparams.dz_res         = ft.length-RPparams.dz_part;
    RPparams.acc_disp_part  = ft.disp*RPparams.dz_part*1e-3;
    RPparams.acc_disp_res   = ft.disp*RPparams.dz_res*1e-3;

    PP(:,k)                 = min_gf(Ein,Eoptfilt,ft,RPparams);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
function gf = min_gf(Ein,Eoptfilt,ft,RPparams)

    P0      = get_power(Ein);

    Ein     = fn(Ein);
    Eoptfilt= fn(Eoptfilt);

    Nepoch  = 1e2;
    gfp     = zeros(Nepoch);
    Costs   = zeros(Nepoch);
    cost_th = 0.1;
    errCost = 1;
    epoch   = 1;

    Costs(1)= get_power(Ein);% mW
    gfp(1)  = ft.gf*1e-6*P0*exp(-ft.alphaLin*RPparams.dz);% 1/(mW.m)

    lr      = 1e-6;
    mr      = 1;

    while Costs(epoch) > cost_th && abs(errCost)> 1e-6 && epoch < Nepoch
%     disp("-----------------------------------")
        [Eprop,Ekerr]   = RPbranch(Ein,ft,RPparams,gfp(epoch));

        Ediff           = diff_fields(Eoptfilt,Eprop);
        Costs(epoch+1)  = get_cost(Eprop,Eoptfilt);
        epoch           = epoch+1;        

        dCost           = get_gradient(Ediff,Ekerr,gfp(epoch-1));
        gfp(epoch)      = update_gf(mr,lr,gfp(epoch-1),dCost);

        errCost         = get_relerr(Costs(epoch-1),Costs(epoch));

    end

    Costs   = Costs(Costs~=0);
    gfp     = gfp(gfp~=0);
    Costs   = Costs(2:end);
    Nepochs = length(Costs)-1;

    subplot(2,2,1)
    scatter(1:Nepochs+1,Costs,"filled")
    title("Cost")
    subplot(2,2,2)
    scatter(1:Nepochs+2,gfp,"filled")
    title("gamma fibre")
    subplot(2,2,[3,4])
    scatter(Costs,gfp(1:end-1),"filled")
    xlabel("cost")
    ylabel("gamma")

    gf = gfp(end);

%     pause(0.5)
%     close all

% ----------------------------------------------------------------------- %
function Cost = get_cost(Eprop,Eoptfilt)
    
    Ediff   = diff_fields(Eprop,Eoptfilt);
    Cost    = get_power(Ediff);

% ----------------------------------------------------------------------- %
function Ediff = diff_fields(E1,E2)
    
    Ediff = E1;
    Ediff.field = E1.field-E2.field;

% ----------------------------------------------------------------------- %
function gf = update_gf(mem_rate,learn_rate,gf,dCost)
    
    old     = mem_rate*gf;
    offset  = -learn_rate*dCost;
    gf      = old+offset;

% ----------------------------------------------------------------------- %
function relerr = get_relerr(varargin)
    
    Xref = varargin{1};
    Xcomp= varargin{2};

    relerr = (Xref-Xcomp)/Xref;
    if nargin >= 3
        disp("absolute relative difference")
        relerr = abs(relerr);
        if nargin == 4
            disp("given in percentage")
            relerr = 100*relerr;
        end
    end

% ----------------------------------------------------------------------- %
function dCost = get_gradient(Ediff,Ekerr,gfp)

    In_Mean = (Ediff.field).*(conj(Ekerr.field));
    In_Real = mean(In_Mean);
    dCost   = -2/gfp*real(In_Real);

    


