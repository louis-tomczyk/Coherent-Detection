function [power_profile,Ecdc_res] = ppe(Ein,Eref,ft,PPEparams)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PPE - POWER PROFILE ESTIMATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-07-07
%   Version         : 1.1
%
% ----- Main idea -----
%   Estimation of the power profile along the fibre using
%   [JLT,TANIMURA,2020]-based algorithm
%
% ----- INPUTS -----

% -------------------------------------------------------------------------
% method.pp = method used to estimate the POWER PROFILE in the fiber.
%           .m = is the MATHEMATICAL operation used. Can be:
%
%               - "cov" : covariance between REF and BP field, in [REF]^2
%
%                       cov(X,Y) = mean([REF-mean(REF)]*.[BP-mean(BP)])
%
%               - "cor" : correlation between REF and BP field, no unit
%
%                                      mean([REF-mean(REF)]*.[BP-mean(BP)])
%                       cor(REF,BP) = ------------------------------------
%                                                   std(REF).STD(BP)
%
%               - "pow" : power in the BP, in [mW]
%                        
%                       pow(BP) = cumtrapz(time,BP)
%
%          .q = is the QUANTITY used for the operation. Can be:
%
%               - "field"   : the normalised fields are used, in [sqrt(mW)]
%               - "power"   : the powers are used, in [mW]
%               - "modulus" : the modulus of the fields, in [sqrt(mW)]
% -------------------------------------------------------------------------
% method.wf = method used for the WAVEFORM reconstruction. Can be:
%           - "tanimura": the estimated non linear phase is:
%
%                                   ΦNL = ε.|u(z,t)|^2
%
%               where: - (u) is fiber output partially CD compensated
%                      - (ε) non linear parameter to empirically tune
% ----------------
%           - "tanimura_losses": the estimated non linear phase is as in
%           "tanimura gamma" plus a power compensation of the field by
%           multiplying with and exponential:
%
%                                   exp(+alphaLin* k*dL/factor)
%
%               where: - (alphaLin) is the linear attenuation, in [1/m]
%                      - (k) is the iterative step number
%                      - (dL) is the elementary length step, in [m]
%                      - (factor) is a factor to be tune empirically to
%                      best match the waveform reconstruction
% ----------------
%           - "true_NLP": the estimated non linear phase is:
%
%                                   ΦNL = gf.|u(z,t)|^2.length_eff_part
%
%               where: - (gf,u) see above
%                      - (length_eff_part) is the partial effective length:
%
%                                                  1-exp(-alphaLin.k.dL)
%                               length_eff_part = ----------------------
%                                                        alphaLin
% ----------------
%           - "true NLP losses": combination of "true NLP" + losses
%           compensation
% -------------------------------------------------------------------------

    power_profile   = zeros(1,PPEparams.link.nsteps_BP);
    Enlc.lambda     = ft.lambda;    

    plotplot = 0;
    if plotplot == 1
        figure
        pause(2)
    end

    Length = PPEparams.link.nsteps_BP*PPEparams.link.dl;

    for k = 1:PPEparams.link.nsteps_BP
        length_part = k*PPEparams.link.dl;
        length_res  = Length-k*PPEparams.link.dl;

        %%% partial CD^{-1}
        Ecdc_part   = cdc(Ein,ft,length_part);
    
        %%% partial NL^{-1}
        if PPEparams.phys.nl_factor ~= 0

            switch PPEparams.method.wf

                case "tanimura"
                    % tanimura way
                    % my way --- fails whatever values of gf from 0.1 ->
                    % 0.0001
%                    PHI_NL = fastexp(mean(0.01*empower(Ecdc_part.field)));
                    % tanimura
                    PHI_NL = fastexp(+PPEparams.phys.nl_factor*empower(Ecdc_part.field));
                    Enlc.field = Ecdc_part.field.*PHI_NL;
%                     fprintf("\n\n")

                case "tanimura_losses"
                    % tanimura way + losses compensation
                    Enlc.field = Ecdc_part.field.*fastexp(+PPEparams.phys.nl_factor*empower(Ecdc_part.field)).*exp(ft.alphaLin*length_part/PPEparams.factor);

                case "NLP_kerr_only"
                    % true non linear phase with effective length
                    length_eff_part = (1-exp(-ft.alphaLin*k*PPEparams.link.dl))/ft.alphaLin;
                    Enlc.field      = Ecdc_part.field.*fastexp(+PPEparams.phys.nl_factor*empower(Eref.field)*length_eff_part);

                case "NLP_kerr_losses"
                    % true non linear phase with effective length + losses compensation
                    length_eff_part = (1-exp(-ft.alphaLin*k*PPEparams.link.dl))/ft.alphaLin;
                    Enlc.field      = Ecdc_part.field.*fastexp(+PPEparams.phys.nl_factor*empower(Eref.field)*length_eff_part).*exp(ft.alphaLin*length_part/PPEparams.factor);
                
                case "NLP_kerr_disp"
                    % true non linear phase with effective length
                    length_eff_part = (1-exp(-ft.alphaLin*k*PPEparams.link.dl))/ft.alphaLin;
                    Enlc.field      = Ecdc_part.field.*fastexp(+PPEparams.phys.nl_factor*empower(Ecdc_part.field)*length_eff_part);

                case "NLP_kerr_disp_losses"
                    % true non linear phase with effective length + losses compensation
                    length_eff_part = (1-exp(-ft.alphaLin*k*PPEparams.link.dl))/ft.alphaLin;
                    Enlc.field      = Ecdc_part.field.*fastexp(+PPEparams.phys.nl_factor*empower(Ecdc_part.field)*length_eff_part).*exp(ft.alphaLin*length_part/PPEparams.factor);
                
                otherwise
                    disp("method not implemented or not recognised")
                    break
            end
 
        else
            Enlc.field = Ecdc_part.field;
        end

        %%% residual CD^{-1}
        Ecdc_res    = cdc(Enlc,ft,length_res);
   

        if plotplot == 1
            clf
            subplot(1,2,1)
                hold on
                scatter(real(Ecdc_part.field),imag(Ecdc_part.field),'MarkerEdgeColor','k',"MarkerFaceColor",'k');
                title(sprintf('%.2f [km]',(PPEparams.link.nsteps_tot-k)*PPEparams.link.dl*1e-3))
                axis([-7,7,-7,7])
            subplot(1,2,2)
                hold on
                scatter(real(Ecdc_res.field),imag(Ecdc_res.field),'MarkerEdgeColor','k',"MarkerFaceColor",'k');
                axis([-7,7,-7,7])
            pause(0.01)
        end

        %%% power profile
        % mathematical estimator selection
        if strcmp(PPEparams.method.pp.m,"cov")==1 || strcmp(PPEparams.method.pp.m,"cor")==1
            if strcmp(PPEparams.method.pp.q,'field') == 1
                X = Ecdc_res.field;
                Y = Eref.field;

            elseif strcmp(PPEparams.method.pp.q,'mod') == 1
                X = abs(Ecdc_res.field);
                Y = abs(Eref.field);

            else
                X = empower(Ecdc_res.field);
                Y = empower(Eref.field);
            end
        else
            X = empower(Ecdc_res.field);
        end

        % quantity on which to apply the estimator
        if strcmp(PPEparams.method.pp.m,"cov")==1

            power_profile(k) = mean(conj(X-mean(X)).*(Y-mean(Y)));
            if strcmp(PPEparams.method.pp.q,'field') == 1
                power_profile(k) = abs(power_profile(k));
            end

        elseif strcmp(PPEparams.method.pp.m,"cor")==1
            tmp = corrcoef(X,Y);
            power_profile(k) = tmp(1,2);
            if strcmp(PPEparams.method.pp.q,'field') == 1
                power_profile(k) = abs(power_profile(k));
            end
        else
            tmp = mean(empower(Ecdc_res.field));
            power_profile(k)    = tmp(end);
        end

    end

end


