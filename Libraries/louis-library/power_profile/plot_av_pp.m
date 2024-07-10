function f = plot_av_pp(PPEparams,plot_params,amp,Mean_PP_trunc,locvals)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : 
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-01-27
%   Version         : 1.1
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    
    f = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        hold on
        plot(PPEparams.plot.dist,Mean_PP_trunc,'k')

        amp_scale   = max(Mean_PP_trunc)-min(Mean_PP_trunc);
        pp_mean     = mean(Mean_PP_trunc);
        ymin        = pp_mean-amp_scale/2;
        ymax        = pp_mean+amp_scale*7/10;
        ylim([ymin,ymax])

        for jj = 0:amp.Nspan
            plot([100,100]*jj,[ymin,ymax],'-.k')
        end

        scatter(locvals.loc_mins*PPEparams.link.dl*1e-3,...
                locvals.val_mins,50,...
                'filled','k')
        scatter(locvals.loc_maxs*PPEparams.link.dl*1e-3,...
                locvals.val_maxs,50,...
                'filled','k')

        xlabel("distance [km]")
        if PPEparams.plot.norm == 0
            ylabel("correlation")
        else
            ylabel("centered and mean-normalised correlation")
        end
        
        Subtitle    = char(plot_params.subtitle);
        Date        = char(datetime);
        year        = Date(8:11);
        index_year  = strfind(Subtitle,year);
        
        % if it starts with 'pd '
        if strcmp(Subtitle(1),'p') == 1
            Subtitle    = Subtitle(1:index_year-11);        
        else
            Subtitle    = Subtitle(7:index_year-11);
        end


        
        title(strcat(plot_params.title," averaged"))
        subtitle(Subtitle)
        
        exportgraphics(f,strcat(plot_params.subtitle,"_averaged.png"),"Resolution",300)
        pause(0.5)
        close all
