function [pp_r,d] = get_ppr(PPEparams,tx,ft,amp,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-03-21
%   Version         : 2.1.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    % compute the ideal power levels
    [~,out_dBm] = prop_PC(tx.PdBm,ft.alphadB,ft.length,2*amp.Nspan,tx.PdBm);

    % set the axis in case of only "light" power map needed
    if strcmp(PPEparams.plot.ref.what,"lightest")==1
        tmp         = zeros(amp.Nspan,2);
    
        for k = 0:amp.Nspan
            tmp(k+1,:) = [k,k+1]*ft.length*1e-3;
        end
        
        Dist_ideal  = tmp;
        Dist_ideal  = reshape(Dist_ideal.',[1,numel(Dist_ideal)]);
    
        for k =3:2:length(Dist_ideal)
            Dist_ideal(k) = Dist_ideal(k)+1;
        end
        
        d       = Dist_ideal(1:end-1);
        pp_r     = out_dBm;
    elseif strcmp(PPEparams.plot.ref.what,"light")==1
        d       = PPEparams.plot.dist;
        nsteps  = floor(PPEparams.link.nsteps_fibre/amp.Nspan);
        pp      = zeros(amp.Nspan,nsteps);

        Pmax    = out_dBm(1);
        Pmin    = out_dBm(2);
        a       = (Pmax-Pmin)/(1-nsteps);
        b       = Pmax -1/(1-nsteps);

        for k = 1:amp.Nspan
            pp(k,:) = a*linspace(1,nsteps,nsteps)+b;
        end

        pp_r = reshape(pp.',1,[]);
    else
        assert(nargin == 6,"AXIS and/or LAS structures as input are missing")
        Axis    = varargin{1};
        las     = varargin{2};
        amp_ref = amp;
        amp_ref = set_topology(tx,ft,amp_ref);
        pp_r    = get_pp(PPEparams,Axis,las,tx,ft,amp_ref);
    end

    if PPEparams.plot.ref.std == 1
        pp_r = standardize(pp_r);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : 
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

function [G_dB,out_dBm] = prop_PC(in_dBm,alpha_dB,L,N,P_wanted_dBm)
    
    out_dBm     = zeros(1,N);
    G_dB        = zeros(1,N);
    out_dBm(1)  = in_dBm;
    
    for k = 2:floor(N/2)+1
        out_dBm(2*(k-1))          = attenuation(out_dBm(2*(k-1)-1),alpha_dB,L);
        [G_dB(k-1),out_dBm(2*k-1)]= EDFA_PC(out_dBm(2*(k-1)),P_wanted_dBm);
    end
    
    G_dB            = G_dB(1:floor(N/2));

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function out_dBm = attenuation(in_dBm,alpha_dB,L,varargin)

    in_W            = dBm2W(in_dBm);
    alpha_lin       = log(10)/10*alpha_dB*1e-3;
    out_W           = in_W*exp(-alpha_lin*L);
    
    if nargin == 4
        ase_level   = varargin{1};
        ASE_W       = dBm2W(ase_level);
        out_W       = out_W + ASE_W;
    end

    out_dBm   = W2dBm(out_W);

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function [G_dB,P_out_dBm] = EDFA_PC(P_in_dBm,P_wanted_dBm)

  P_in_W      = dBm2W(P_in_dBm);
  P_wanted_W  = dBm2W(P_wanted_dBm);
  P_ase_dBm   = ase(P_in_dBm,50);
  P_ase_W     = dBm2W(P_ase_dBm);
  
  G_lin       = (P_wanted_W-P_ase_W)./P_in_W;
  G_dB        = 10*log10(G_lin);
  P_out_W     = G_lin*P_in_W + P_ase_W;
  P_out_dBm   = W2dBm(P_out_W);
  
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function out_dBm = ase(in_dBm,offset_dB)

    out_dBm = in_dBm-offset_dB;

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------    
function out = dBm2W(in)

    out = 10.^(in/10)*1e-3;

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPR - Power Profile Reference
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.com
%   Date            : 2023-01-27
%   Version         : 2.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function out = W2dBm(in)

    out = 10*log10(in/1e-3);


