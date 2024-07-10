function export_pp(PP,plot_params)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : EXPORT_PP - EXPORT POWER PROFILE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-13
%   Version         : 1.2.1
%
% ----- MAIN IDEA -----
%   Export the power profiles.
%
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    if plot_params.norm == 1
        pathh = strcat(pwd,'/',plot_params.alpha,'/txt/stdized/',plot_params.what_charac,'/pp');
    else
        pathh = strcat(pwd,'/',plot_params.alpha,'/txt/not stdized/',plot_params.what_charac,'/pp');
    end
    
    cd(pathh)

    PP_raw_pd       = transpose(PP{1}{1});
    PP_mean_pd      = transpose(PP{2}{1});
    PP_no_pd        = transpose(PP{3}{1});
    if size(PP,2) == 4
        PP_mean_no_pd   = transpose(PP{4}{1});
    end

    PP_raw_pd_c     = clean_data(PP_raw_pd,'row',4);
    PP_mean_pd_c    = clean_data(PP_mean_pd,'row',4);
    PP_no_pd_c      = clean_data(PP_no_pd,'row',4);

    if size(PP,2) == 4
        PP_mean_no_pd_c = clean_data(PP_mean_no_pd,'row',4);
    end

    writematrix(PP_raw_pd,...
    strcat(plot_params.subtitle,' --- PP_raw_pd ','.csv'));

    writematrix(PP_mean_pd,...
    strcat(plot_params.subtitle,' --- PP_av_pd ','.csv'));

    writematrix(PP_no_pd,...
    strcat(plot_params.subtitle,' --- PP_av_no_pd ','.csv'));
    
    if size(PP,2) == 4
        writematrix(PP_mean_no_pd,...
        strcat(plot_params.subtitle,' --- PP_moved_av_no_pd ','.csv'));
    end

    %%% 
    writematrix(PP_raw_pd_c,...
    strcat(plot_params.subtitle,' --- PP_raw_pd_cleaned','.csv'));

    writematrix(PP_mean_pd_c,...
    strcat(plot_params.subtitle,' --- PP_av_pd_cleaned','.csv'));

    writematrix(PP_no_pd_c,...
    strcat(plot_params.subtitle,' --- PP_av_no_pd_cleaned','.csv'));
    
    if size(PP,2) == 4
        writematrix(PP_mean_no_pd_c,...
        strcat(plot_params.subtitle,' --- PP_moved_av_no_pd_cleaned','.csv'));
    end
    
    cd ~/Documents/6___Telecom_Paris/3_Codes/louis/Optilux/Data/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_out = clean_data(data_in,row_or_col,step)

    % values for ROW_OR_COL
    %   - ROW --- keep one row every 2 rows
    %   - COL --- keep one col every 2 cols
    
    % values for ODD_OR_EVEN
    %   - ODD --- start from 1st ROW/COL
    %   - EVEN --- start from 2nd ROW/COL
    
    if strcmp(row_or_col,'row') == 1
        data_out = data_in(2:step:end,:);
    elseif strcmp(row_or_col,'col') == 1
        data_out = data_in(1,2:step:end,:);
    end
    

