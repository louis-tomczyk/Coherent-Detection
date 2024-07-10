function [Eout,amp] = channel(Ein,ft,amp)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CHANNEL
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@enst.fr
%   Date            : 2023-01-27
%   Version         : 2.0
%
% ----- Main idea -----
%   Simulate consecutive lossy spans
%
% ----- INPUTS -----
%   EIN     (sructure)  The field to propagate
%   FT      (structure) The fibre parameters
%   AMP     (structure) The amplifiers parameters
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    n_spans = length(amp.losses.ids);

    for k=1:n_spans
        mm          = amp.losses.ids(k);
        losses_tmp  = amp.losses.(strcat('span',num2str(mm)));
        dloss       = losses_tmp(:,1);
        loss_dB     = losses_tmp(:,2);
    
        if k == 1
%             disp("span No 1")
            Eout = choose_amp_mode(Ein,ft,dloss,loss_dB,amp);
%             get_power(Eout,struct('unit','dBm','polar','tot'))
        else
%             sprintf("span No %i",k)
            Eout = choose_amp_mode(Eout,ft,dloss,loss_dB,amp);
%             get_power(Eout,struct('unit','dBm','polar','tot'))
        end
    end

    % update the number of spans for potential following identical
    % ideal fibres
%     amp.Nspan = amp.Nspan - n_spans;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CHOOSE_AMP_MODE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@enst.fr
%   Date            : 2023-01-18
%   Version         : 1.0
%
% ----- Main idea -----
%   Select which kind of lossy span is wanted: with PC or GC amp mode
%
% ----- INPUTS -----
%   EIN     (sructure)  The field to propagate
%   FT      (structure) The fibre parameters
%   DLOSS   (array)[km] The locations of the losses
%   LOSS_DB (array)[dB] The value of the loss
%   AMP     (structure) The amplifiers parameters
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

function Eoutfib = choose_amp_mode(Ein,ft,dloss,loss_dB,amp)

    if strcmp(amp.mode,'GC') == 1
        assert(strcmp(amp.mode,'GC')==1,"should be Gain Control mode")
        Eoutfib = lossy_span(Ein,ft,dloss,loss_dB,amp);
    else
        assert(strcmp(amp.mode,'PC')==1,"should be Power Control mode")
        PC_dBm  = amp.outpower;%get_power(Ein,struct('unit',"dBm",'polar','tot'));
        Eoutfib = lossy_span(Ein,ft,dloss,loss_dB,amp,PC_dBm);
    end


% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : LOSSY SPAN - see SPAN_WITH_LOSS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@enst.fr
%   Date            : 2023-01-17
%   Version         : 1.1
%
% ----- Main idea -----
%   Create a SNGLE lossy span
%
% ----- INPUTS -----
%   EIN     (sructure)  The field to propagate
%   FT      (structure) The fibre parameters
%   DLOSS   (array)[km] The locations of the losses
%   LOSS_DB (array)[dB] The value of the loss
%   AMP     (structure) The amplifiers parameters
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
%   Functions           : FIBER - AMPLIFLAT
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------
function [Eout,amp] = lossy_span(Ein,ft,dloss,loss_dB,amp,varargin)


%     ptmp = get_power(Ein,struct('unit','dBm','polar','tot'));
%     sprintf("power before [fibre+loss] = %.1f",ptmp)

    assert(length(dloss) == length(loss_dB),...
        ['the number of distance of losses should ' ...
        'match the number of loss values'])

    if nargin == 6
        % check that we are in Power Constant mode
        assert(strcmp(amp.mode,'PC')==1)
        if strcmp(amp.mode,'PC') ==1
            tmp             = amp;
            tmp.outpower    = varargin{1};
        end
    else
        % check that we are in Gain Constant mode
        assert(strcmp(amp.mode,'GC')==1, "The expected power is missing")
        tmp = amp;
    end

    %%% propagation in the lossy fibre before (included) the last loss
    ft_tmp  = ft;

    if loss_dB ~= 0
        for k = 1:length(dloss)
            if k == 1
                ft_tmp.length    = dloss(1);
            else
                if k <= length(dloss)
                    ft_tmp.length    = dloss(k)-dloss(k-1);
                else
                    ft_tmp.length = ft.length-dloss(end);
                end
            end
    %         ft_tmp.length
            loss_lin    = log(10)/10*loss_dB(k)/2;
    
            %%% aplying the loss to the fibre
            if k == 1
                Eout    = fiber(Ein,ft_tmp);
            else
                Eout    = fiber(Eout,ft_tmp);
            end
    
    %         ptmp = get_power(Eout,struct('unit','dBm','polar','tot'));
    %         sprintf("power before [loss] = %.1f",ptmp)
            
            Eout.field  = Eout.field*exp(-loss_lin);
    %         ptmp = get_power(Eout,struct('unit','dBm','polar','tot'));
    %         sprintf("power after [loss] = %.1f",ptmp)
        end
    
    
        %%% propagation into the last part of the span
        ft_tmp.length   = ft.length-dloss(end);
        Eout            = fiber(Eout,ft_tmp);
    %     ptmp            = get_power(Eout,struct('unit','dBm','polar','tot'));
    %     sprintf("power at end of the span = %.1f",ptmp)
    
        %%% amplification at the end of the span
        if strcmp(amp.mode,'PC') == 1
            if isfield(tmp,'gain') == 1
                tmp = rmfield(tmp,"gain");
            end
        end
    else
        Eout        = fiber(Ein,ft);
    end

    if strcmp(amp.mode,'PC') == 1
        if isfield(tmp,'gain') == 1
            tmp = rmfield(tmp,"gain");
        end
    end
    Eout          = ampliflat(Eout,tmp);
%     ptmp = get_power(Eout,struct('unit','dBm','polar','tot'))
%     sprintf("power at end of the span's amplifier = %.1f",ptmp)


        