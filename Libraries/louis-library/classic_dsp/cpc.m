function Ecpc = cpc(Ein,tx,rx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CPC - CARRIER PHASE COMPENSATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-05
%   Version         : 1.3
%
% ----- Main idea -----
%   Carrier phase (frequency + phase noise) compensation using
%   VITERBI-VITERBI algorithm.
%   This code comes, 99.99% of it, from OPTILUX 2009 - DSP4COHDEC
%   I only take a small of part of it to put it in a new function.
%
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to be compensated
%           structure should be organised as:
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
%   NAVG:   (scalar)[] number of points used to averge the phase noise
%   MODFOR: (string) modulation format
%
% ----- OUTPUTS -----
%  ECPEC:   (structure) containing the Carrier Phase Estimation
%           Compensated field.
%           structure containing
%               - LAMBDA (scalar)[nm]: wavelength
%               - FIELD (array)[sqrt(mW)]: normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions   : DSP4COHDEC
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2009
%   Title of program    : Optilux
%   Code version        : 2009
%   Type                : Optical simulator toolbox - source code
%   Web Address         : Partage Zimbra "Phd Louis Tomczyk"/Optilux/
%                           Optilux v2009
% ---------------------------------------------

    modfor  = tx.modfor;
    navg    = rx.CPC.navg;
    
    % valence number of the modulation format
    switch modfor
        case "bpsk"
            M   = 2;
        case "qpsk"
            M   = 4;
        case "dqpsk"
            M   = 4;
        case "16qam"
            M   = 16;
        otherwise
            disp("not implented yet")
    end

    if navg ~= 0
        % Estimating frequency in bell labs style [citation needed]:
        omega = cumsum(vitvit(Ein.field.*conj(fastshift(Ein.field,1)),M,M,navg,false));

%         figure
%             subplot(1,3,1)
%             hold all
%             plot(unwrap(angle(Ein.field(:,1))))
%             plot(unwrap(angle(Ein.field(:,2))))
%             legend("X","Y")
%             title("louis")
%     
%             subplot(1,3,2)
%             hold all
%             plot(omega(:,1))
%             plot(omega(:,2))
%             legend("X","Y")
%             title(sprintf("optilux before - navg = %i",navg))

        % Cleaning omega to match the circularity:
        closestallowedendpoints = omega(1,:)+round((omega(end,:)-omega(1,:))/2/pi)*2*pi;
        correctionratio         = closestallowedendpoints./omega(end,:);
    
        omega = ((omega-ones(length(omega),1)*omega(1,:)).*(ones(length(omega),1)*correctionratio))...
            +ones(length(omega),1)*omega(1,:);

%             subplot(1,3,3)
%             hold all
%             plot(omega(:,1))
%             plot(omega(:,2))
%             legend("X","Y")
%             title("optilux after")
%         legend("myway","before","after")

        navg    = 3;
        theta   = vitvit(Ein.field,P,M,navg,true);
    
        if P > 1
            CarrierPhaseOffSet = +pi/4;
        else
            CarrierPhaseOffSet = 0;
        end
    
        Carrier = fastexp(-theta+CarrierPhaseOffSet);
    end
    
    Ecpc.lambda= Ein.lambda;
    Ecpc.field = Ein.field.*Carrier;
    Ecpc.Nsps = length(Ecpc.field)/tx.Nsymb;
    Ecpc = fn(Ecpc);
    
    if isfield(tx,'Plin')
        Ecpc.field = Ecpc.field*sqrt(tx.Plin/2);
    end

return
        % Demodulating signals:
        sigdemod = Ein.field.*fastexp(-omega);
   
        % Estimating phase using Viterbi and Viterbi method:
        navg    = 5;
        P       = 2;
        theta   = vitvit(sigdemod,P,M,navg,true);
    
        if P > 1
            CarrierPhaseOffSet = +pi/4;
        else
            CarrierPhaseOffSet = 0;
        end
    
        Carrier = fastexp(-omega-theta+CarrierPhaseOffSet);
    
    else
        P       = 2;
        navg    = 3;
        theta   = vitvit(Ein.field,P,M,navg,true);
    
        if P > 1
            CarrierPhaseOffSet = +pi/4;
        else
            CarrierPhaseOffSet = 0;
        end
    
        Carrier = fastexp(-theta+CarrierPhaseOffSet);
    end
    
    Ecpc.lambda= Ein.lambda;
    Ecpc.field = Ein.field.*Carrier;
    
    Ecpc = fn(Ecpc);
    
    if isfield(tx,'Plin')
        Ecpc.field = Ecpc.field*sqrt(tx.Plin/2);
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = vitvit(s,P,M,k,applyunwrap)
    L = length(s);
    % Phase -> Phase * M
    % Amplitude -> Amplitude ^ P

    if P == M
        s = s .^ P;
    else
        s = abs(s).^P .* fastexp( angle( s.^ M ) );
    end

    if k>0
        N = 2*k + 1;
        if N<L
            size(ones(1,size(s,2)));
            size(fft( ones(N, 1) / N , L ));
            Smoothing_Filter    = fft( ones(N, 1) / N , L )*ones(1,size(s,2));
            % usefull for TRUE DATA with ENSSAT data formating. Examples with QPSK1.mat
            % in TP CD ENSSAT folder
%             Smoothing_Filter    = fft( ones(N, 1) / N , L )'*ones(1,size(s,2))';
            s                   = ifft( fft(s,L) .* Smoothing_Filter );
        else
            slong               = repmat( s, ceil(N / L), 1 );
            Smoothing_Filter    = fft( ones(N, 1) / N , ceil(N / L).*L ) * ones(1,size(s,2));
            slong               = ifft( fft(slong) .* Smoothing_Filter );
            s                   = slong(1:L,:);
        end
    end

    if applyunwrap
        theta = unwrap( angle(s) ) / M;

    else
        theta = angle( s ) / M;
    end
return