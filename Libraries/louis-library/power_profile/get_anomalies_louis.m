function [pps,losshat] = get_anomalies(ppref,pploss,amp,PPEparams,ft,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_ANOMALIES
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-31
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Estimate the losses from power profile estimations
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    %%% Profiles scaling for surimposition
    Mean_lossy          = mean(pploss);
    Mean_ref            = mean(ppref);
    
    diff_means          = Mean_ref-Mean_lossy;
    
    if diff_means>0
        pploss_scaled   = pploss-diff_means;
    else
        pploss_scaled   = pploss+diff_means;
    end
    
    %%% Anomalies - localisation estimation
    ppdiff              = abs(pploss_scaled-ppref);
    z                   = PPEparams.plot.dist;
        
    [z_1, pploss_1]     = get_derivative(z,pploss,1);
    pploss_1            = standardize(pploss_1);
    [phat,zhat]         = max(ppdiff);
    zhat                = PPEparams.plot.dist(zhat);

    res                 = PPEparams.link.dl/1000;
    resolution          = get_decimals(res);
    zhat                = round(zhat,resolution.n);

    %%% Anomalies - extrema estimation
    tmpMx               = zeros(amp.Nspan,1);
    tmpMy               = zeros(amp.Nspan,1);
    tmpmx               = zeros(amp.Nspan,1);
    tmpmy               = zeros(amp.Nspan,1);

    dl                  = PPEparams.link.dl;
    ratio_dl            = dl/1000;
    N_steps_span        = PPEparams.link.nsteps_fibre/amp.Nspan;

    % this IF loop to check that the extrema estimation is working
    % correctly
    if nargin == 6 && varargin{1} == 1
%         figure('Color',[0 0 0])
%         hold on
%         plot(z_1,pploss_1,'k','LineWidth',2)
    
        for k = 1:amp.Nspan
            if k == amp.Nspan
%                 plot(z_1(1+(k-1)*N_steps_span:k*N_steps_span-1), ...
%                     pploss_1(1+(k-1)*N_steps_span:k*N_steps_span-1));

                [tmpmy(k),tmpmx(k)] = min(pploss_1(1+(k-1)*N_steps_span-1:k*N_steps_span-1));
                tmpmx(k)            = (k-1)*ft.length/dl+tmpmx(k);
            else
%                 plot(z_1(1+(k-1)*N_steps_span:k*N_steps_span), ...
%                     pploss_1(1+(k-1)*N_steps_span:k*N_steps_span));

                [tmpmy(k),tmpmx(k)] = min(pploss_1(1+(k-1)*N_steps_span:k*N_steps_span));
                tmpmx(k)            = (k-1)*ft.length/dl+tmpmx(k);
                [tmpMy(k),tmpMx(k)] = max(pploss_1(1+(k-1)*N_steps_span:k*N_steps_span));
                tmpMx(k)            = (k-1)*ft.length/dl+tmpMx(k);
            end
    
%             scatter(tmpmx(k)*ratio_dl,tmpmy(k),100,'filled')
%             scatter(tmpMx(k)*ratio_dl,tmpMy(k),100,'filled')
   
        end
    else   
        for k = 1:amp.Nspan
            if k == amp.Nspan
                [tmpmy(k),tmpmx(k)] = min(pploss_1(1+(k-1)*N_steps_span-1:k*N_steps_span-1));
                tmpmx(k)            = (k-1)*ft.length/dl+tmpmx(k);
            else
                [tmpmy(k),tmpmx(k)] = min(pploss_1(1+(k-1)*N_steps_span:k*N_steps_span));
                tmpmx(k)            = (k-1)*ft.length/dl+tmpmx(k);
                [tmpMy(k),tmpMx(k)] = max(pploss_1(1+(k-1)*N_steps_span:k*N_steps_span));
                tmpMx(k)            = (k-1)*ft.length/dl+tmpMx(k);
            end  
        end
    end

    %%% Outputs
    pps.loss.der0       = pploss;
    pps.loss.der0_scaled= pploss_scaled;
    pps.loss.der1       = pploss_1;
    pps.loss.diff       = ppdiff;

    losshat.loc         = zhat;
    losshat.val.mins.x  = tmpmx;
    losshat.val.mins.y  = tmpmy;
    losshat.val.maxs.x  = tmpMx;
    losshat.val.maxs.y  = tmpMy;
    losshat.val.amps    = losshat.val.maxs.y(1:end-1)-losshat.val.mins.y(2:end);
    

    %%% Plot the outputs
    if nargin == 6 && varargin{1} == 1
    
        figure('units', 'normalized', 'outerposition', [0 0 1 1],'Color',[1,1,1]);
            subplot(1,3,1)
                hold all
                plot(z,ppref,'k','LineWidth',3)
                plot(z,pploss,'color','c','LineWidth',3)
                plot(z,pploss_scaled,'b','LineWidth',1)
                xlabel(sprintf("distance [km] - resolution = %.1f [km]",res))
                ylabel("standardized estimated power [a.u.]")
                legend("ref",'monitored','scaled monitored')
                title("Power profile estimation")
            subplot(1,3,2)
                hold on
                plot(z,ppdiff,'k','LineWidth',2)
                scatter(zhat,phat,100,"MarkerFaceColor",'b',"MarkerEdgeColor",'b')
    
                text(zhat+25,phat*1.02,sprintf("zhat = %.1f [km]",zhat))
                xlabel(sprintf("distance [km] - resolution = %.1f [km]",res))
                legend("|scaled lossy- reference|",'location','northwest')
                title("Anomaly detection - localisation of loss")
            subplot(1,3,3)
                hold on
                plot(z,pploss,'b','LineWidth',2)
                plot(z_1,pploss_1,'k','LineWidth',2)
    
                scatter(tmpmx*ratio_dl,tmpmy,100,'filled')
                scatter(tmpMx*ratio_dl,tmpMy,100,'filled')
                xlim([10,780])
                ylim([min(pploss_1),max(pploss_1)]*1.1)
                xlabel(sprintf("distance [km] - resolution = %.1f [km]",res))
                legend("pploss",'diff pploss')
                title("Anomaly detection - amplitude of loss")

    end