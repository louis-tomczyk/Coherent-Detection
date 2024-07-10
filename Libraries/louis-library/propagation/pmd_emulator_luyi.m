function [Eout,ft] = pmd_emulator_luyi(Ein,Axis,ft)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PMD_EMULATOR
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 
%   Date            : 2024-01-29: creation
%   Version         : 1.0
%
% ----- Main idea -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
%   Functions           : PMD_EMULATOR
%   Author              : Elie AWWAD
%   Author contact      : elie.awwad@telecom-paris.fr
%   Date                : 2023-03-21
%   Title of program    : NA
%   Code version        : NA
%   Type                : Phenix
%   Web Address         : NA
% ---------------------------------------------

% ft      = set_channel_params(ft,Axis);
% ft      = plates_definition(ft);
% ft      = create_channel(ft);
% Eout    = apply_pmd(Ein,ft);

[V,h11,h12,h21,h22,Y] = pmd_emulator(Ein.field,Axis.freq,ft.pmdpar,0,1);

Eout.field = Y;
ft.Hplates = V;
ft.Hmimo.HH = h11;
ft.Hmimo.HV = h12;
ft.Hmimo.VH = h21;
ft.Hmimo.VV = h22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------
% ----- CONTENTS -----
%   ft = set_channel_params(ft,Axis)
%   ft = create_channel(ft)
%   [Eout,rx]   = equalisation_simu(E,Axis,rx)
%   plot_constellation_in(E,iter,isample,rx)
% ---------------------------------------------

function ft = set_channel_params(ft,Axis)

    if isfield(ft,"pdl") == 0
        ft.pdl = 0;
    end

    ft.tauPMD   = ft.pmdpar/sqrt(ft.nplates);
    ft.gamma    = 10^(-0.05*ft.pdl/sqrt(ft.nplates));
    ft.exp_tau  = exp(1i*pi.*Axis.freq.*ft.tauPMD);
    
    ft.Hmimo0.HH = ones(size(Axis.freq));
    ft.Hmimo0.HV = zeros(size(Axis.freq));
    ft.Hmimo0.VH = ft.Hmimo0.HV;
    ft.Hmimo0.VV = ft.Hmimo0.HH;
% ---------------------------------------------

function ft = plates_definition(ft)

    if isfield(ft,"polMatrix") == 0

        % 6 = 3 (PMD) + 3(PDL)
        % PMD / PDL == Rout*H*Rin
%         for k = 1:6:6*ft.nplates
%             V{:,:,k}    = zeros(2);
%             V{:,:,k+2}  = zeros(2);
%             V{:,:,k+3}  = zeros(2);
%             V{:,:,k+5}  = zeros(2);
%         end

        for k = 1:6:6*ft.nplates
            theta       = 2*pi*rand;% %2*pi*rand;
            V{1,1,k}    = cos(theta);
            V{1,2,k}    = sin(theta);
            V{2,1,k}    = -V{1,2,k};
            V{2,2,k}    = V{1,1,k};
    
		    psi = 2*pi*rand;				
            V{1,1,k+1} = ft.exp_tau*exp(1i*psi/2);
            V{1,2,k+1} = 0;
            V{2,1,k+1} = 0;
            V{2,2,k+1} = conj(ft.exp_tau)*exp(-1i*psi/2);
    
            V{1,1,k+2} = V{1,1,k};
            V{1,2,k+2} = V{2,1,k};
            V{2,1,k+2} = V{1,2,k};
            V{2,2,k+2} = V{2,2,k};
    
            theta = 2*pi*rand;
            V{1,1,k+3} = cos(theta);
            V{1,2,k+3} = sin(theta);
            V{2,1,k+3} = -V{1,2,k+3};
            V{2,2,k+3} = V{1,1,k+3};
    
            V{1,1,k+4} = 1;
            V{1,2,k+4} = 0;
            V{2,1,k+4} = 0;
            V{2,2,k+4} = ft.gamma;
    
            V{1,1,k+5} = V{1,1,k+3};
            V{1,2,k+5} = V{2,1,k+3};
            V{2,1,k+5} = V{1,2,k+3};
            V{2,2,k+5} = V{2,2,k+3};
        end
    

    else
        return
    end

    ft.Hmimo = V;
% ---------------------------------------------

function ft = create_channel(ft)
    
    % 1-3   = 7-9   = 13-15 = ... PMD Rotation
    %  2    =  8    =   14  = ... PMD Delay

    % 4-6   = 10-12 = 16-18 = ... PDL Rotation
    %  5    =   11  =   17  = ... PDL Delay

    h11 = ft.Hmimo0.HH;
    h12 = ft.Hmimo0.HV;
    h21 = ft.Hmimo0.VH;
    h22 = ft.Hmimo0.VV;

    for k = 1:6*ft.nplates
    
        O11 = ft.Hmimo{1,1,k};
        O12 = ft.Hmimo{1,2,k};
        O21 = ft.Hmimo{2,1,k};
        O22 = ft.Hmimo{2,2,k};
    
        X11 = O11.*h11 + O12.*h21;
        X12 = O11.*h12 + O12.*h22;
        X21 = O21.*h11 + O22.*h21;
        X22 = O21.*h12 + O22.*h22;
    
        h11 = X11;
        h12 = X12;
        h21 = X21;
        h22 = X22;

    end

    ft_tmp = ft;
    ft_tmp = rmfield(ft_tmp,'Hmimo');
    clear ft
    ft = ft_tmp;

    ft.Hmimo.HH = h11;
    ft.Hmimo.HV = h12;
    ft.Hmimo.VH = h21;
    ft.Hmimo.VV = h22;
% ---------------------------------------------

function Eout = apply_pmd(Ein,ft)

    XHft              = FFT(Ein.field(:,1)).';
    XVft              = FFT(Ein.field(:,2)).';

    YHft            = ft.Hmimo.HH.*XHft+ft.Hmimo.HV.*XVft;
    YVft            = ft.Hmimo.VH.*XHft+ft.Hmimo.VV.*XVft;

    Eout.field(:,1) = iFFT(YHft);
    Eout.field(:,2) = iFFT(YHft);
    Eout.lambda     = Ein.lambda;
% ---------------------------------------------

% =========================== %
% ---- SUB SUB FUNCTIONS ---- %
% =========================== %

function [Rin,Rout] = get_rotation_matrix()

    theta   = 2*pi*rand;
    Rin     = [cos(theta), sin(theta);-sin(theta), cos(theta)];
    Rout    = Rin(:,:,1)';
% ---------------------------------------------










%% TRASH
% psi                     = 2*pi*rand;	
% [V{:,:,k},V{:,:,k+2}]   = get_rotation_matrix();
% [V{:,:,k+3},V{:,:,k+5}] = get_rotation_matrix();
% 
% V{:,:,k+1}              = [ft.exp_tau*exp(+1i*psi/2),0;
%                            0,conj(ft.exp_tau)*exp(-1i*psi/2)];
% V{:,:,k+4}              = [1,0; 0,ft.gamma];














% 
% tmpHH   = ft.Hmimo{1,1,k};
% tmpHV   = ft.Hmimo{1,2,k};
% tmpVH   = ft.Hmimo{2,1,k};
% tmpVV   = ft.Hmimo{2,2,k};
% 
% if k == 1
%     XHH = tmpHH.*ft.Hmimo0.HH + tmpHV.*ft.Hmimo0.VH;
%     XHV = tmpHH.*ft.Hmimo0.HV + tmpHV.*ft.Hmimo0.VV;
%     XVH = tmpVH.*ft.Hmimo0.HH + tmpVV.*ft.Hmimo0.VH;
%     XVV = tmpVH.*ft.Hmimo0.HV + tmpVV.*ft.Hmimo0.VV;
% else
%     XHH = tmpHH.*ft.Hmimo.HH + tmpHV.*ft.Hmimo.VH;
%     XHV = tmpHH.*ft.Hmimo.HV + tmpHV.*ft.Hmimo.VV;
%     XVH = tmpVH.*ft.Hmimo.HH + tmpVV.*ft.Hmimo.VH;
%     XVV = tmpVH.*ft.Hmimo.HV + tmpVV.*ft.Hmimo.VV;
% end
% 
% ft.HmimoTot.HH   = XHH;
% ft.HmimoTot.HV   = XHV;
% ft.HmimoTot.VH   = XVH;
% ft.HmimoTot.VV   = XVV;