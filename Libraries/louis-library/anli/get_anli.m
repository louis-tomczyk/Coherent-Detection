function [Enli,nli_stats] = get_anli(Ein,Ecpc,Axis,tx,rx,ft)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_ANLI
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-12
%   Version         : 1.4
%
% ----- Main idea -----
%   Get the Non Linear Interference (NLI) coefficient statistical
%   properties following the article procedure.
%   They should have mentionned that the DEMODULATION operation is
%   DATA-AIDED which means, it uses the phases values of each symbol.
%   
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to be normlised
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
% ----- BIBLIOGRAPHY -----
%   Functions           : Opt_filt & Elec_filt
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : Unknown
%   Title of program    : Code DP_QPSKand16QAM
%   Code version        : 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% -----------------------
%   Articles
%   Author              : Nicolas ROSSI, Petros RAMANTANIS, Jean-Claude
%                         ANTONA
%   Title               : Nonlinear Interference Noise Statistics in
%                         Unmanagged Coherent Networks with Channnels
%                         Propagating over Differennt Lightpaths
%   Jounal              : ECOC
%   Volume - N°         : Cannes
%   Date                : 2014
%   DOI                 : 10.1109/ECOC.2014.6964043
% ---------------------------------------------

    % CPC has Nsymb-Ntaps (from CMA) symbols
    Nsymb_cpc           = size(Ecpc.field(:,1));

    % optical and electrical filtering
    Einprep             = rx_filt(Ein,Axis,tx);

    % compensating the pre dispersion if there is
    Rx                  = rx;
    Rx.CDC.len_tot      = tx.pd/ft.disp*1e3;
    Einprep_cdc         = cdc(Einprep,ft,Rx);

    % downsample to 1SPS and match the number of symbols
    Einds               = ds(Einprep_cdc,1,Axis.Nsps,1);
    Einds.field(Nsymb_cpc+1:end,:)  = [];

    % normalisation (ISAIA is not normalising)
    Eindsnorm           = fn(Einds,tx.Plin);
    Ecpcnorm            = fn(Ecpc,tx.Plin);

    %%% demodulation - data aided version
    % it will bring all of the points with null imaginary part, (no phase)
    [Eindemod,Phi_in]   = qpsk_demod(Eindsnorm,Axis);
    [Ecpcdemod,Phi_cpc] = qpsk_demod(Ecpcnorm,Axis);
    Ecpcdemod = qpsk_demod_distance(Ecpcnorm,Axis)
    [Phicorr,shifts]    = correct_phase(Phi_in,Phi_cpc,Axis,rx.CMA);

    % then remodulate with the corrected phase PHI, see CORRECT_PHASE
    % function's help
    Ecpcmod.lambda      = Ecpcdemod.lambda;
    Ecpcmod.field(:,1)  = Ecpcdemod.field(:,1).*fastexp(Phicorr.pol_1);
    Ecpcmod.field(:,2)  = Ecpcdemod.field(:,2).*fastexp(Phicorr.pol_2);    

    % We can then re-demodulate using the input phase
    Ecpcdemod           = qpsk_demod(Ecpcmod,Phi_in);

    Ecpcdemod           = qpsk_demod_distance(Ecpcnorm,Axis);
    % to finally truncate the fields to remove [PHI_CPC_a ... PHI_CPC_b],
    % see CORRECT_PHASE function's help
    Ecpctrunc                       = Ecpcdemod;
    Ecpctrunc.field(1:shifts.x,1)   = NaN;
    Ecpctrunc.field(1:shifts.y,2)   = NaN;
    tmp_cpc                         = rmmissing(Ecpctrunc.field);
    Ecpcfinal.lambda                = Ecpctrunc.lambda;
    Ecpcfinal.field                 = tmp_cpc;

    % in case it fails, the point is centered at (0,1), then we turn it to
    % the right place
    mean_angle_X = mean(angle(Ecpcfinal.field(:,1)));
    mean_angle_Y = mean(angle(Ecpcfinal.field(:,2)));
    
    if mean_angle_X > 1
        Ecpcfinal.field(:,1)= (Ecpcfinal.field(:,1)).*fastexp(-pi/2); % +pi/4 if polarisation coupling
    end
    if mean_angle_Y > 1
        Ecpcfinal.field(:,2)= (Ecpcfinal.field(:,2)).*fastexp(-pi/2); % +pi/4 if polarisation coupling    
    end

    Eindemod.field(1:shifts.x,1)    = NaN;
    Eindemod.field(1:shifts.y,2)    = NaN;
    tmp_in                          = rmmissing(Eindemod.field);
    Einfinal.lambda                 = Ecpctrunc.lambda;
    Einfinal.field                  = tmp_in;

    % get the field of noise
    % ROSSI way
    Enli        = diff_fields(Ecpcfinal,Einfinal);

%     % ISAIA way
%     Enli = Ecpcfinal;
%     Enli.field(:,1) = Enli.field(:,1)-mean(Enli.field(:,1));
%     Enli.field(:,2) = Enli.field(:,2)-mean(Enli.field(:,2));

    % get the statistics
    nli_stats   = get_stats(Enli,"var","b");
    
    if strcmp(rx.NLN.stats,'all') == 1
        nli_stats.skw   = get_stats(Enli,'skw','b');
        nli_stats.kts   = get_stats(Enli,'kts','b');
    end

    % plot if wanted
    if rx.NLN.plot == 1
        Ein     = Eindsnorm;
        Ecpc    = Ecpcnorm;
        plot_const(Ein,Einfinal,Ecpc,Ecpcfinal);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------
function [phi,shifts] = correct_phase(Phi_cpc,Phi_in,Axis,CMAparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CORRECT_PHASE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-15
%   Version         : 1.4
%
% ----- Main idea -----
%   Correct the phase values after CPC for constellation demodulation.
%   Let's write PHI_IN the input phase used for modulate the signal of N
%   samples:
%
%       PHI_IN  = [ PHI_IN_1 ... PHI_IN_c | PHI_IN_Shift ... PHI_IN_N]
%
%   If the Carrier Compensation (CC) is effective then the phases should
%   almost match (between INPUT and CC) with the original values and the
%   sequence will be shifted because of the moving average.
%   The first values will be then PHI_IN_1, PHI_IN_2, upto PHI_IN_c (the
%   'c' letter does not stand for anything. Same for the letters 'a' and 
%   'b' later on.
%
%   The objective is to find how much shifted it was.
%   Let's write PHI_CPC the phase after the CC:
%
%       PHI_CPC = [ PHI_IN_Shift ... PHI_IN_N | PHI_CPC_a ... PHI_CPC_b]
%
%   Let's write PHI the corrected phase
%
%       PHI     = [ PHI_CPC_a ... PHI_CPC_b | PHI_IN_Shift ... PHI_IN_N]
%
% ----- INPUTS -----
% ----- BIBLIOGRAPHY -----
%   Functions           :
%   Author              :
%   Author contact      :
%   Date                :
%   Title of program    :
%   Code version        :
%   Type                :
%   Web Address         :
% -----------------------
%   Articles
%   Author              :
%   Title               :
%   Jounal              :
%   Volume - N°         :
%   Date                :
%   DOI                 :
% ---------------------------------------------

    n_polar = length(fieldnames(Phi_in));

    if n_polar == 1
        XCORRx      = abs(xcorr(Phi_in.pol_1,Phi_cpc.pol_1));
        assert(sum(XCORRx==fliplr(XCORRx)) == size(XCORRx,1),"The cross-correlation should be perfectly symetric." + ...
            "PMDC not good enough.")

        [~,ix]      = max(XCORRx);
        [~,ix_flip] = max(fliplr(XCORRx));
        assert(ix==ix_flip,"The cross-correlation should be perfectly symetric." + ...
            "maxima do not match.")

        diff_x      = abs(ix-Axis.Nsymb);
        shifts.x    = diff_x-CMAparams.taps;
        phi.pol_1   = Phi_in.pol_1(diff_ix+1:end-CMAparams.taps+dix);
        phi.pol_1   = circshift(Phi_in.pol_1,shifts.x);
    else
        % the cross correlation is used to find where the signals overlap
        % the most
        XCORRx      = abs(xcorr(Phi_in.pol_1,Phi_cpc.pol_1));
        XCORRy      = abs(xcorr(Phi_in.pol_2,Phi_cpc.pol_2));

        % assert that the classical DSP steps are working well
        assert(sum(XCORRx==fliplr(XCORRx)) == size(XCORRx,1),"The cross-correlation should be perfectly symetric." + ...
            "PMDC or CPC not good enough on the X-polar")
        assert(sum(XCORRy==fliplr(XCORRy)) == size(XCORRy,1),"The cross-correlation should be perfectly symetric." + ...
            "PMDC or CPC not good enough on the Y-polar")

        % maxima location of the cross-correlation
        [~,ix]      = max(XCORRx);
        [~,iy]      = max(XCORRy);
        
        % assert that flipping the cross correlation does the same maxima
        % location. Should be the case if perfectly symetric.
        [~,ix_flip] = max(fliplr(XCORRx));
        [~,iy_flip] = max(fliplr(XCORRy));

        assert(ix==ix_flip,"The cross-correlation should be perfectly symetric." + ...
            "correlation maximum do not match on X-polar.")
        assert(iy==iy_flip,"The cross-correlation should be perfectly symetric." + ...
            "correlation maximum do not match on Y-polar.")

        % Location in the sequence of where the shift (PHI_IN_Shift) is.
        diff_x      = abs(ix-Axis.Nsymb);
        diff_y      = abs(iy-Axis.Nsymb);
        
        % phase shift due to CPC
        shifts.x    = diff_x-CMAparams.taps;
        shifts.y    = diff_y-CMAparams.taps;

        % shifting the phases
        phi.pol_1   = circshift(Phi_in.pol_1,shifts.x);
        phi.pol_2   = circshift(Phi_in.pol_2,shifts.y);

    end
                    
%-----------------------------------------------------
function Ediff = diff_fields(E1,E2)

    Ediff       = E1;
    Ediff.field = E1.field-E2.field;
