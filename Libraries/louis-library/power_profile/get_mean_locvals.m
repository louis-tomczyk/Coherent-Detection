function locvals_m = get_mean_locvals(locvals,what)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_MEAN_LOCVALS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-15
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   The power profile estimator script MAIN_PP writes in output files the
%   characteristic values of the profile obtained.
%   This function enables to GET the (x) LOCATIONS and (y) VALUES according
%   to whether the span number or the parameter under study values
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

    %% maintenance

    loc_mins = locvals.loc_mins;
    loc_maxs = locvals.loc_maxs;

    val_mins = locvals.val_mins;
    val_maxs = locvals.val_maxs;
    
    %% calculations
    if strcmp(what,"perspan") == 1
        %%% locs
        locvals_m.perspan.loc_mins.mean   = transpose(mean(loc_mins.'));
        locvals_m.perspan.loc_mins.std    = transpose(std(loc_mins.'));
        
        locvals_m.perspan.loc_maxs.mean   = transpose(mean(loc_maxs.'));
        locvals_m.perspan.loc_maxs.std    = transpose(std(loc_maxs.'));
        
        %%% values
        locvals_m.perspan.val_mins.mean   = transpose(mean(val_mins.'));
        locvals_m.perspan.val_mins.std    = transpose(std(val_mins.'));
        
        locvals_m.perspan.val_maxs.mean   = transpose(mean(val_maxs.'));
        locvals_m.perspan.val_maxs.std    = transpose(std(val_maxs.'));
        
    elseif strcmp(what,"percharac") == 1
        %%% locs
        locvals_m.percharac.loc_mins.mean = transpose(mean(loc_mins));
        locvals_m.percharac.loc_mins.std  = transpose(std(loc_mins));
        
        locvals_m.percharac.loc_maxs.mean = transpose(mean(loc_maxs));
        locvals_m.percharac.loc_maxs.std  = transpose(std(loc_maxs));
        
        %%% values
        locvals_m.percharac.val_mins.mean = transpose(mean(val_mins));
        locvals_m.percharac.val_mins.std  = transpose(std(val_mins));
        
        locvals_m.percharac.val_maxs.mean = transpose(mean(val_maxs));
        locvals_m.percharac.val_maxs.std  = transpose(std(val_maxs));

    elseif strcmp(what,"both") == 1

        locvals_m.perspan   = get_mean_locvals(locvals,"perspan");
        locvals_m.percharac = get_mean_locvals(locvals,"percharac");

        tmp_perspan         = locvals_m.perspan.perspan;
        locvals_m           = rmfield(locvals_m,"perspan");
        locvals_m.perspan   = tmp_perspan;

        tmp_percharac       = locvals_m.percharac.percharac;
        locvals_m           = rmfield(locvals_m,"percharac");
        locvals_m.percharac = tmp_percharac;

    end
end