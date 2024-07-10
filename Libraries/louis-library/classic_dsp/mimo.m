function [Emimo,rx] = mimo(Ecdc,Axis,rx)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : MIMO - MULTI INPUT MULTI OUTPUT
%                    previously: PMDC - POLARISATION MODE DISPERSION COMPENSATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   ArXivs          : 2024-01-11: cleaning Yves' code
%                   : 2024-01-17: including Élie's code
%                   : 2024-01-18: adding Lingliu singularity management algorithm
%   Date            : 2024-01-29: possibility to change learning rate for equalisation
%   Version         : 2.1.1
%
% ----- Main idea -----
%   Polarisation Mode Dispersion compensation using CMA algorithm.
%
% ----- INPUTS -----
%   ECDC    (structure) of the Chromatic Dispersion Compensated field
%               - LAMBDA    (scalar)[nm]        Wavelength
%               - FIELD     (array)[sqrt(mW)]   Normalised electric fields
%   Axis    (structure) axis parameters, see SET_AXIS
%   RX      (structure) recever parameters, see SET_RX
%
% ----- OUTPUTS -----
%  EPMDC    (structure) of the Chromatic Dispersion Compensated field
%               - LAMBDA    (scalar)[nm]        Wavelength
%               - FIELD     (array)[sqrt(mW)]   Normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions           : SOPRECOVERY_FSE_CMAANDDRE
%   Author              : Yves JAOUEN
%   Author contact      : simu.jaouen@telecom-paris.fr
%   Date                : 2021-02-09
%   Title of program    : NA
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------
%   Functions           : Polarization Demultiplexer
%   Author              : Élie AWWAD
%   Author contact      : phenix.awwad@telecom-paris.fr
%   Date                : 2023-06-19
%   Title of program    : Phenix
%   Code version        : 2023/10
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    Ein = fn(Ecdc);

    if strcmp(rx.CMA.method,"simu") == 1

        [Etrain,rx]     = init_simu(Ein,Axis,rx);
        rx              = train_simu(Etrain,rx);
        [Emimo,rx]      = equalisation_simu(Etrain,Axis,rx);

    elseif strcmp(rx.CMA.method,"phenix") == 1

        [Etrain,rx]     = init_phenix(Ein,rx);
        rx              = train_phenix(Etrain,rx);
        [Emimo,rx]      = equalisation_phenix(Ein,rx);
    end

    rx = is_mimo_ok(Emimo,rx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED FUNCTIONS - YVES JAOUEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
% --- from YVES JAOUEN
%   [E,rx]      = init_simu(Ecdc,Axis,rx)
%   rx          = train_simu(E,rx)
%   [Eout,rx]   = equalisation_simu(E,Axis,rx)
%   plot_constellation_in(E,iter,isample,rx)
% 
% --- from ÉLIE AWWAD
%   [Etrain,Eq,rx]  = init_phenix(E,rx)
%   E               = Default_MIMO_Equalizer_in(rx)
%   [Eq,rx]         = train_phenix(Etrain,Eq,rx)
%   rx              = check_singularity(Eq,rx)
%   [Emimo,rx]      = equalisation_phenix(E,Eq,rx)
%   d               = downsampleOTS_in( x, r, delay)
%   plot_fir_const(Emimo,rx)
%   rx              = extract_loss(Ein,rx)
% ---------------------------------------------

function [Etrain,rx] = init_simu(Ein,Axis,rx)

    if rx.CMA.Ntrain == Axis.Nsamp
        powerOf2        = nextpow2(Axis.Nsymb);
        rx.CMA.Ntrain   = 2^(powerOf2-1);
    end

    Etrain               = (Ein.field).';

    if isfield(Axis,"ADC") == 0
        Axis.ADC.Nsamp  = 2*Axis.Nsymb;
    end

    Nsps_in         = Axis.ADC.Nsamp/Axis.Nsymb;
    assert(abs(round(Nsps_in,4)- rx.CMA.Nsps_in)<1e-3,'---  there should be 2 samples/symbol at this step.')

    halfTap         = floor(rx.CMA.taps/2)+1;
    rx.CMA.Nplot    = 1000;

    if isfield(rx.CMA,"Ntrain") == 0
        if Axis.Nsamp > 20000
            rx.CMA.Ntrain = ceil((Axis.ADC.Nsamp-rx.CMA.taps)/2)-1;
        else
            rx.CMA.Ntrain = floor(Axis.Nsamp/2) - rx.CMA.Nplot - ceil(rx.CMA.Nsps_in);
        end
    end

    rx.CMA.FIR.HH            = zeros(1,rx.CMA.taps);
    rx.CMA.FIR.HH(1,halfTap) = 1.0 ;
    rx.CMA.FIR.VV            = rx.CMA.FIR.HH;
    rx.CMA.FIR.HV            = zeros(1,rx.CMA.taps);
    rx.CMA.FIR.VH            = rx.CMA.FIR.HV;
% ---------------------------------------------

function rx = train_simu(E,rx)

    rx.CMA.FIR.lossH = zeros(1,rx.CMA.Ntrain);
    rx.CMA.FIR.lossV = zeros(1,rx.CMA.Ntrain);

    for k = 1:rx.CMA.Ntrain

         istart = 1+(k-1)*rx.CMA.Nsps_in;
         iend   = istart+rx.CMA.taps-1;

         Hin    = E(1,istart:iend);
         Vin    = E(2,istart:iend);

         Xout   = rx.CMA.FIR.HH*Hin.' + rx.CMA.FIR.HV*Vin.';
         Yout   = rx.CMA.FIR.VH*Hin.' + rx.CMA.FIR.VV*Vin.';
     
         rx.CMA.FIR.lossH(k)    = 2*(abs(Xout)^2 - 1)*Xout;
         rx.CMA.FIR.lossV(k)    = 2*(abs(Yout)^2 - 1)*Yout;

         rx.CMA.FIR.HH          = rx.CMA.FIR.HH - rx.CMA.lr_train*rx.CMA.FIR.lossH(k)*conj(Hin);
         rx.CMA.FIR.HV          = rx.CMA.FIR.HV - rx.CMA.lr_train*rx.CMA.FIR.lossH(k)*conj(Vin);
         rx.CMA.FIR.VH          = rx.CMA.FIR.VH - rx.CMA.lr_train*rx.CMA.FIR.lossV(k)*conj(Hin);
         rx.CMA.FIR.VV          = rx.CMA.FIR.VV - rx.CMA.lr_train*rx.CMA.FIR.lossV(k)*conj(Vin);

        if rx.CMA.plot == 1
            plot_constellation_in(E,iter,isample,rx)
        end
    end 

    if rx.CMA.plot == 1
        figure(Name="CMA cost function")
            subplot(1,2,1)
            loglog(smoothdata(abs(rx.CMA.FIR.lossH),'gaussian', 100))
            ylim([1e-3,10])
            subplot(1,2,2)
            loglog(smoothdata(abs(rx.CMA.FIR.lossV),'gaussian', 100))
            ylim([1e-3,10])
    end
% ---------------------------------------------

function [Eout,rx] = equalisation_simu(E,Axis,rx)
    
    if isfield(Axis,"ADC") == 0
        Axis.ADC.Nsamp  = 2*Axis.Nsymb;
    end
    
    rx.CMA.Nsymbols = floor(Axis.ADC.Nsamp/rx.CMA.Nsps_in) - rx.CMA.taps;
    Eout            = zeros(2, Axis.ADC.Nsamp/rx.CMA.Nsps_in);
    
    for k = 1:rx.CMA.Nsymbols
         
         istart     = 1+(k-1)*rx.CMA.Nsps_in;
         iend       = istart+rx.CMA.taps-1;
         
         Hin        = E(1,istart:iend);
         Vin        = E(2,istart:iend);

         Eout(1,k)  = rx.CMA.FIR.HH*Hin.' + rx.CMA.FIR.HV*Vin.';
         Eout(2,k)  = rx.CMA.FIR.VH*Hin.' + rx.CMA.FIR.VV*Vin.'; 
    end  

    tmp = Eout;
    clear Eout
    Eout.lambda = 1550;
    Eout.field  = tmp.';
% ---------------------------------------------

function plot_constellation_in(E,iter,isample,rx)

    x = iter/rx.CMA.Ntrain;
     try
         if rx.CMA.plot == 1
             % plotting the reak-time evolution of the constellation
             % using 1000 points every 500 iterations
             if (mod(iter,rx.CMA.Nplot) == 1)
                istart = isample;
    
                Eout = zeros(2,rx.CMA.Nplot);

                for k = 1:rx.CMA.Nplot
                    iend        = istart+rx.CMA.taps-1;

                    Hin         = E(1,istart:iend);
                    Vin         = E(2,istart:iend);  
    
                    Eout(1,k)   = rx.CMA.FIR.HH*Hin.' + rx.CMA.FIR.HV*Vin.';
                    Eout(2,k)   = rx.CMA.FIR.VH*Hin.' + rx.CMA.FIR.VV*Vin.';
    
                    istart      = istart+rx.CMA.Nsps_in;
                end 
                 
                figure(2)
                
                    subplot(2,2,1)
                        polarplot(angle(Eout(1,:)), abs(Eout(1,:)),'.')
                        title(sprintf('After %3d symbols', iter'))

                    subplot(2,2,2)
                        polarplot(angle(Eout(2,:)), abs(Eout(2,:)),'.')
                        title(sprintf('After %3d symbols', iter'))

                    subplot(2,2,3)
                        hold on
                        scatter(1:rx.CMA.taps,abs(rx.CMA.FIR.HH),...
                            "MarkerEdgeColor",[0,0.5,1]*(1-4*x),"MarkerFaceColor",[0,0.5,1]*(1-4*x));
                        scatter(1:rx.CMA.taps,-abs(rx.CMA.FIR.HV),...
                            "MarkerEdgeColor",[1,1,1]*(1-2*x),"MarkerFaceColor",[1,1,1]*(1-2*x));

                        title('hxx and hxy')
                        axis([1,rx.CMA.taps,-1,1])

                    subplot(2,2,4)
                        hold on
                        scatter(1:rx.CMA.taps,abs(rx.CMA.FIR.VV),...
                            "MarkerEdgeColor",[0,0.5,1]*(1-4*x),"MarkerFaceColor",[0,0.5,1]*(1-4*x));
                        scatter(1:rx.CMA.taps,-abs(rx.CMA.FIR.VH),...
                            "MarkerEdgeColor",[1,1,1]*(1-2*x),"MarkerFaceColor",[1,1,1]*(1-2*x));

                        xlabel('Taps number')
                        title('hyy and hyx')
                        axis([1,rx.CMA.taps,-1,1])
                        
                pause(0.1)
             end
        end
     catch
    
     end
% ---------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED FUNCTIONS - ÉLIE AWWAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Etrain,rx] = init_phenix(E,rx)
    
    % FIR init, output  2D
    rx          = Default_MIMO_Equalizer_in(rx);

    % training set
    if rx.CMA.Ntrain > size(E.field,2)
        rx.CMA.Ntrain = size(E.field,1);
    end

    Etrain.field = E.field(1:rx.CMA.Ntrain,:);
% ---------------------------------------------

function rx = train_phenix(Etrain,rx)

    %%% 1st try
    [Etmp.field,rx.CMA.FIR]  = poldemux_core(Etrain.field,rx.CMA.FIR, ...
                            rx.CMA.lr_train,rx.CMA.taps,rx.CMA.Rv,rx.CMA.R,rx.CMA.step,0);
    

    %%% singularity management
    rx = check_singularity(rx);

    if rx.CMA.singularity.flag == true
        flag_count  = 1;

        while flag_count < rx.CMA.singularity.loops
            
            %%% plot status
            plot_fir_const(Etmp,rx);
            if rx.CMA.plot == 1
                set(gcf, 'name', sprintf("training loop = %d",flag_count))
            end

            %%% re-init the taps
            rx = Default_MIMO_Equalizer_in(rx);
            
            %%% update the taps
            [Etmp.field,rx.CMA.FIR]    = poldemux_core(Etrain.field,rx.CMA.FIR, ...
                                   rx.CMA.lr_train,rx.CMA.taps,rx.CMA.Rv,rx.CMA.R,rx.CMA.step,0);

            %%% re-check singularity
            rx = check_singularity(rx);

            %%% if singularity, re-train
            flag_count = flag_count+1;
            if rx.CMA.singularity.flag == false
                 break
            end

        end
    end

    Etmp = ds(Etmp,2,1,1);
    plot_fir_const(Etmp,rx)

    if rx.CMA.plot_final == 1
        set(gcf, 'name', "CMA training: end")
    end

    rx = reshapeFIR2D(rx);
% ---------------------------------------------

function [Emimo,rx] = equalisation_phenix(E,rx)

    [Emimo.field,rx.CMA.FIR] = poldemux_core(E.field,rx.CMA.FIR, ...
        rx.CMA.lr_eq, rx.CMA.taps, rx.CMA.Rv, rx.CMA.R, rx.CMA.step, 0);

    Emimo           = ds(Emimo,2,1,1);

    plot_fir_const(Emimo,rx)
    if rx.CMA.plot_final == 1
        set(gcf, 'name', "CMA equalised")
    end
% ---------------------------------------------

function plot_fir_const(Emimo,rx)

    if rx.CMA.plot_final == 1 || rx.CMA.singularity.plot == 1

        TAPS    = (1:rx.CMA.taps)-ceil(rx.CMA.taps/2);
        t       = linspace(-1,1,100);
        yc1     = sqrt(1-t.^2);
        yc2     = -sqrt(1-t.^2);

        figure;
            subplot(2,2,1)
                hold on
                plot(t,yc1,'--',t,yc2,'--')
                scatter(real(Emimo.field(:,1)),imag(Emimo.field(:,1)), ...
                        10,MarkerEdgeColor="k",MarkerFaceColor="k")
                axis([-1,1,-1,1]*sqrt(2))
                pbaspect([1 1 1])
                title("Pol H",FontSize=15,FontWeight="bold")
                set(gca,"FontSize",15)
            subplot(2,2,2)
                hold on
                plot(t,yc1,'--',t,yc2,'--')
                scatter(real(Emimo.field(:,2)),imag(Emimo.field(:,2)), ...
                        10,MarkerEdgeColor="k",MarkerFaceColor="k")
                axis([-1,1,-1,1]*sqrt(2))

                pbaspect([1 1 1])
                title("Pol V",FontSize=15,FontWeight="bold")
                set(gca,"FontSize",15)

            if isfield(rx.CMA.FIR,"HH") == 0
                tmp             = rx.CMA.FIR;
                rx.CMA          = rmfield(rx.CMA,'FIR');

                if length(size(tmp)) ~= 3
                    tmp = reshape(tmp,[rx.CMA.taps,2,2]);
                end

                rx.CMA.FIR.HH    = tmp(1:rx.CMA.taps,1,1);
                rx.CMA.FIR.HV    = tmp(1:rx.CMA.taps,2,1);
                rx.CMA.FIR.VH    = tmp(1:rx.CMA.taps,1,2);
                rx.CMA.FIR.VV    = tmp(1:rx.CMA.taps,2,2);

                clear tmp
            end
            subplot(2,2,3)
                hold on
                plot(TAPS,abs(rx.CMA.FIR.HH),'k',LineWidth=2)
                plot(TAPS,abs(rx.CMA.FIR.HV),'--k',LineWidth=1.5)
                ylim([0,1]*sqrt(2))
                legend("$HH$","$HV$",FontWeight="bold",FontSize=13)
                set(gca,"FontSize",15)
            subplot(2,2,4)
                hold on
                plot(TAPS,abs(rx.CMA.FIR.VH),'--k',LineWidth=1.5)
                plot(TAPS,abs(rx.CMA.FIR.VV),'k',LineWidth=2)
                ylim([0,1]*sqrt(2))
                legend("$HV$","$VV$",FontWeight="bold",FontSize=13)
                set(gca,"FontSize",15)
    end
% ---------------------------------------------

% =========================== %
% ---- SUB SUB FUNCTIONS ---- %
% =========================== %


function rx = Default_MIMO_Equalizer_in(rx)

    if rx.CMA.singularity.flag == true

        if length(size(rx.CMA.FIR)) == 2
            rx = reshapeFIR3D(rx);
        end

        if strcmp(rx.CMA.singularity.method,"rand") == 1
            psi = pi*rand();
            azi = pi*rand();
    
            % mat = 
            %  a  b
            % -b* a*
            mat = ...
                 [cos(psi)*cos(azi)-1i*sin(psi)*sin(azi), ...
                 -sin(psi)*cos(azi)+1i*cos(psi)*sin(azi); ...
                 sin(psi)*cos(azi)+1i*cos(psi)*sin(azi), ...
                 cos(psi)*cos(azi)+1i*sin(psi)*sin(azi) ]; 
            rx.CMA.FIR = zeros(rx.CMA.taps, 2, 2);
            rx.CMA.FIR(round(rx.CMA.taps/2),:,:) = mat;

        elseif strcmp(rx.CMA.singularity.method,"lingliu") == 1

            FIR        = rx.CMA.FIR;
            HH         = FIR(:,1,1);
            VH         = FIR(:,1,2);
            
            FIR(:,2,1) = -conj(fliplr(VH));
            FIR(:,2,2) = conj(fliplr(HH));

            rx.CMA.FIR = FIR;
        end
    else
        if strcmp(rx.CMA.method_init,"dirac") == 1
            rx.CMA.FIR = zeros(rx.CMA.taps, 2, 2);
            rx.CMA.FIR(round(rx.CMA.taps/2),:,:) = eye(2);
        elseif strcmp(rx.CMA.method_init,"rand") == 1
            psi = pi*rand();
            azi = pi*rand();
            mat = ...
                 [cos(psi)*cos(azi)-1i*sin(psi)*sin(azi), ...
                 -sin(psi)*cos(azi)+1i*cos(psi)*sin(azi); ...
                 sin(psi)*cos(azi)+1i*cos(psi)*sin(azi), ...
                 cos(psi)*cos(azi)+1i*sin(psi)*sin(azi) ]; 
            rx.CMA.FIR = zeros(rx.CMA.taps, 2, 2);
            rx.CMA.FIR(round(rx.CMA.taps/2),:,:) = mat;
        else
            disp("MIMO > INIT_PHENIX: please choose between {dirac,rand} options")
            return
        end
    end

    if length(size(rx.CMA.FIR)) == 3
        rx = reshapeFIR2D(rx);
    end
% ---------------------------------------------

function rx = check_singularity(rx)

    rx      = reshapeFIR3D(rx);
    sum_E   = sum(rx.CMA.FIR,1);
    Det     = abs(det(squeeze(sum_E)));

    if Det < rx.CMA.singularity.thresh
        disp("CMA singularity")
        rx.CMA.singularity.flag = true;
        rx.CMA.singularity.H    = squeeze(sum_E);
        rx.CMA.singularity.Det  = abs(Det);
    else
        if rx.CMA.singularity.flag == true
            disp("CMA singularity - disappeared")
        end
        rx.CMA.singularity.flag = false;
    end

    rx = reshapeFIR2D(rx);
% ---------------------------------------------

function rx = reshapeFIR2D(rx)

    if length(size(rx.CMA.FIR)) == 3
        rx.CMA.FIR = reshape(flip(rx.CMA.FIR,1),[numel(rx.CMA.FIR),1]);
    end

    rx.CMA = sort_struct_alphabet(rx.CMA);
% ---------------------------------------------

function rx = reshapeFIR3D(rx)

    if length(size(rx.CMA.FIR)) == 2
        rx.CMA.FIR = reshape(rx.CMA.FIR,[rx.CMA.taps,2,2]);
    end

    rx.CMA = sort_struct_alphabet(rx.CMA);
% ---------------------------------------------

function rx = is_mimo_ok(Emimo,rx)
% perfect compensation should give equal power 
% in both polarisation and equal to 1/2

    % reseting the flags to avoid miscounting
    if isfield(rx.CMA.flag,'polH')
        rx.CMA.flag = rmfield(rx.CMA.flag,'polH');
    end
    if isfield(rx.CMA.flag,'polV')
        rx.CMA.flag = rmfield(rx.CMA.flag,'polV');
    end
    if isfield(rx.CMA.flag,'count_fail') == 0
        rx.CMA.flag.count_fail = 0;
    end
    if isfield(rx.CMA.flag,'count_singularity') == 0
        rx.CMA.flag.count_singularity = 0;
    end

    % error calculations
    Pmimo   = get_power(Emimo);
    tmpH    = round(Pmimo(2)-1,4);
    tmpV    = round(Pmimo(3)-1,4);

    % concatenating errors if several draws
    if isfield(rx.CMA.flag,"errH") == 0 ||isfield(rx.CMA.flag,"errV") == 0
        rx.CMA.flag.errH = tmpH;
        rx.CMA.flag.errV = tmpV;
    else
        rx.CMA.flag.errH = [rx.CMA.flag.errH,tmpH];
        rx.CMA.flag.errV = [rx.CMA.flag.errV,tmpV];
    end 

    % flag raising
    if abs(tmpH) > rx.CMA.flag.thresh && abs(tmpV) > rx.CMA.flag.thresh
        rx.CMA.flag.polH = "fail";
        rx.CMA.flag.polV = "fail";

    elseif abs(tmpH) > rx.CMA.flag.thresh && abs(tmpV) < rx.CMA.flag.thresh
        rx.CMA.flag.polH = "fail";

    elseif abs(tmpV) > rx.CMA.flag.thresh && abs(tmpH) < rx.CMA.flag.thresh
        rx.CMA.flag.polV = "fail";
    end

    if isfield(rx.CMA.flag,"polH") == 1 || isfield(rx.CMA.flag,"polV") == 1
        disp(rx.CMA.flag)
    end

    % counting the failures
    if isfield(rx.CMA.flag,"polH") == 1 || isfield(rx.CMA.flag,"polV") == 1
        rx.CMA.flag.count_fail = rx.CMA.flag.count_fail+1;
    end

    if rx.CMA.singularity.flag == 1
        rx.CMA.flag.count_singularity = rx.CMA.flag.count_singularity+1;
    end
% ---------------------------------------------




%% TRASH
% function rx = extract_loss_in(Ein,rx)
% 
%     Nsamp = length(Ein.field);
%     
%     assert(nargin == 2, "Missing the Receiver's structure")
%     if isfield(rx.CMA,'loss') == 0
%         rx.CMA.loss.Navg    = 100;
%         rx.CMA.loss.method  = "bloc";
%     end
% 
%     rx.CMA.loss.Nsteps  = floor(Nsamp/rx.CMA.loss.Navg);
%     rx.CMA.loss.LossH   = zeros(rx.CMA.loss.Nsteps,1);
%     rx.CMA.loss.LossV   = zeros(rx.CMA.loss.Nsteps,1);
% 
%     for k = 1:rx.CMA.loss.Nsteps
%         if strcmp(rx.CMA.loss.method,"bloc") == 1
%             tmpEh = Ein.field(1+(k-1)*rx.CMA.loss.Navg:k*rx.CMA.loss.Navg,1);
%             tmpEv = Ein.field(1+(k-1)*rx.CMA.loss.Navg:k*rx.CMA.loss.Navg,2);
%         elseif strcmp(rx.CMA.loss.method,"moving avg") == 1
%             tmpEh = Ein.field(k:k+rx.CMA.loss.Navg,1);
%             tmpEv = Ein.field(k:k+rx.CMA.loss.Navg,2);
%         end
% 
%         tmpEh = abs(tmpEh).^2;
%         tmpEv = abs(tmpEv).^2;
% 
%         rx.CMA.loss.polH(k) = mean((tmpEh-1).^2);
%         rx.CMA.loss.polV(k) = mean((tmpEv-1).^2);
%     end
% 
%     rx.CMA.loss.tot = rx.CMA.loss.polH+rx.CMA.loss.polV;
% % ---------------------------------------------


% function d = downsampleOTS_in( x, r, delay)
% % Function Y = DOWNSAMPLEOTS(X, R, DELAY)
% %   Downsamples an input vector X taking one every R samples. If DELAY is
% %   given it starts from the Y(1) = X(DELAY), Y(2)=X(DELAY+R),....
% %   By default DELAY = 1;
% 
% if nargin < 3   %si le nombre d'argument de la fonction est < 3
%     delay = 1;
% end
% 
% if r == 1
%     d = x(delay:end,:);
%     return
% end
% 
% % Il faut s'assurer que la taille de x est un multiple de r
% over_length = mod(length(x),r);
% x           = x(1:end-over_length,:);
% 
% if length(delay) > 1
%     d1 = x( delay(1) : r : end , 1);
%     d2 = x( delay(2) : r : end , 2);
%     % Si le pic est different pour les 2 polars, il se peut que l'on
%     % recupere un echantillon de plus sur une d'elles. Il faut alors le
%     % supprimer pour egaliser les longueurs de sequences.
%     if length(d1) > length(d2)
%         d(:,1) = d1(1:end-1,1);
%         d(:,2) = d2(:,1);
%     elseif length(d2) > length(d1)
%         d(:,1) = d1(:,1);
%         d(:,2) = d2(1:end-1,1);
%     else
%         d(:,1) = d1(:,1);
%         d(:,2) = d2(:,1);
%     end
% else
%     d = x( delay : r : end , :);
% end
% % ---------------------------------------------
