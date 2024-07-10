function [y,loc] = flmm(Mean_PP,Mean_PPav,what,PPEparams,amp)

    % FIND LOCAL MINIMA AND MAXIMA
    N_points_max = ceil(PPEparams.findloc.uncertainty/PPEparams.link.dl);

    %% minima location
    % --- first location
    if strcmp(what,"min")

        % 1/ find all local minima in both PP
        loc_av  = find(islocalmin(Mean_PPav));
        N_loc_av= length(loc_av);

        x = zeros(1,length(loc_av));
        y = zeros(1,length(loc_av));
    
        loc = zeros(1,N_loc_av);
        
        for k = 1:length(loc_av)
            if loc_av(k)+N_points_max > length(Mean_PP)
                N_points_max_sup = length(Mean_PP)-loc_av(k);
                X = loc_av(k)-N_points_max:loc_av(k)+N_points_max_sup;
            else
                if loc_av(k)+N_points_max+5 > length(Mean_PP)
                    X = loc_av(k)-N_points_max+5:length(Mean_PP);
                else
                    X = loc_av(k)-N_points_max+5:loc_av(k)+N_points_max+5;
                end
            end

            [y(k),x(k)] = min(Mean_PP(X));
            loc(k) = loc_av(k)+x(k)-floor(N_points_max/2+2);

            % last point in the profile should be the end of the last span
            if loc(k) > length(Mean_PP)
                loc(k) = length(Mean_PP);
            end

        end
        
    %% maxima location
    % --- first location
    else
        loc_av  = find(islocalmax(Mean_PPav));
        N_loc_av= length(loc_av);
       
        % If no bell shape curve in the first span (strictly decreasing)
        % then the maximum is the first element
        tmp         = Mean_PPav;
        tmp_diff    = diff(tmp);

        if tmp_diff(1)<0
            loc = zeros(1,N_loc_av+1);
            x   = zeros(1,length(loc_av)+1);
            y   = zeros(1,length(loc_av)+1);
        else
            loc = zeros(1,N_loc_av);
            x   = zeros(1,length(loc_av));
            y   = zeros(1,length(loc_av));
        end

        for k = 1:length(loc_av)
            % start of the link
            if loc_av(k)-N_points_max < 1
                X = 1:loc_av(k)+N_points_max;
            else
                X = loc_av(k)-N_points_max+5:loc_av(k)+N_points_max+5;
            end

            [y(k),x(k)] = max(Mean_PP(X));
            loc(k)      = loc_av(k)+x(k)-floor(N_points_max/2+2);
            y(k)        = Mean_PP(loc(k));

            if loc(k) < PPEparams.link.nsteps_fibre/amp.Nspan
                loc(k) = loc(k)+4;
            end
        end

        if tmp_diff(1)<0
            loc(end)    = 1;
            y(end)      = Mean_PP(1);

            loc         = circshift(loc,1);
            y           = circshift(y,1);
        end

    end
end
