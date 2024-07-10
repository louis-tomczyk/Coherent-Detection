function  [f,t] = plot_raw_pp(PPEparams,plot_params,PP_all,varargin)

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
    if PPEparams.repet.Ntries > 1
        % --- back propagation only in the link
        if strcmp(PPEparams.link.method_BP,'no_pd') == 1

            Mean_PP             = varargin{1};
            Mean_PP_av          = varargin{2};

            hold on
            plot(PPEparams.plot.dist,PP_all,'k')
            plot(PPEparams.plot.dist,Mean_PP)
            plot(PPEparams.plot.dist,Mean_PP_av)
           
        else % --- back propagation in link + predisp
            Mean_PP_trunc       = varargin{1};
            Mean_PP_trunc_av    = varargin{2};

            t = tiledlayout(2,1,'TileSpacing','Compact');
            nexttile,
                plot(PP_all,'k')
            nexttile,
                hold on
                plot(PPEparams.plot.dist,Mean_PP_trunc,'--b')
                plot(PPEparams.plot.dist,Mean_PP_trunc_av,'r')

                amp_scale   = max(Mean_PP_trunc)-min(Mean_PP_trunc);
                pp_mean     = mean(Mean_PP_trunc);
                ymin        = pp_mean-amp_scale/2;
                ymax        = pp_mean+amp_scale*4/3;
                ylim([ymin,ymax])

            legend("truncated PP","truncated averaged PP",'location',"northwest")
        end
    else % Ntries = 1
        if strcmp(PPEparams.link.method_BP,'no_pd') == 1

            PP_av               = varargin{1};

            hold on
            plot(PPEparams.plot.dist,PP_all,'k')
            plot(PPEparams.plot.dist,PP_av)
    
        else % --- back propagation in link + predisp

            PP_trunc            = varargin{1};
            PP_trunc_av         = varargin{2};

            t = tiledlayout(2,1,'TileSpacing','Compact');
            nexttile,
                plot(PP_all,'k')
            nexttile,
                hold on
                plot(PPEparams.plot.dist,PP_trunc)
                plot(PPEparams.plot.dist,PP_trunc_av)

                amp_scale   = max(PP_trunc)-min(PP_trunc);
                pp_mean     = mean(PP_trunc);
                ymin        = pp_mean-amp_scale/2;
                ymax        = pp_mean+amp_scale*4/3;
                ylim([ymin,ymax])
                
    
        end
    end
    
    xlabel(t,"distance [km]")
    if PPEparams.plot.norm == 1
        ylabel(t,"centered and mean-normalised correlation")
    else
        ylabel(t,"correlation")
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

    title(t,strcat(plot_params.title," raw"))
    subtitle(t,Subtitle)

    exportgraphics(f,strcat(plot_params.subtitle,"_raw.png"),"Resolution",300)
    pause(0.5)
%     close all
end