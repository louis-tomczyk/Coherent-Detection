function slopes = get_slopes(locvals,what)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_SLOPES
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-14
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get the characteristic slopes in the profile from characteristic
%   locations and values, whether per span or per parameter under study
%   values.
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
    %%% maintenance

    %%% maintenance
    loc_mins            = locvals.loc_mins;
    loc_maxs            = locvals.loc_maxs;
    val_mins            = locvals.val_mins;
    val_maxs            = locvals.val_maxs;
    
    %%% calculations
    if strcmp(what,"negative")
        dvals = val_mins-val_maxs;
        dlocs = loc_mins-loc_maxs;

        dvals = replace_zeros(dvals);

        Slopes = dvals./dlocs;
        slopes.(what).percharac.mean   = mean(Slopes).';
        slopes.(what).percharac.std    = std(Slopes).';
        slopes.(what).perspan.mean     = mean(transpose(Slopes)).';
        slopes.(what).perspan.std      = std(transpose(Slopes)).';
        
    elseif strcmp(what,"positive")
        val_max = val_maxs(2:end,:);
        val_min = val_mins(1:end-1,:);
        loc_max = loc_maxs(2:end,:);
        loc_min = loc_mins(1:end-1,:);

        dvals   = val_max-val_min;
        dlocs   = loc_max-loc_min;

        dvals = replace_zeros(dvals);

        Slopes = dvals./dlocs;
        slopes.(what).percharac.mean   = mean(Slopes).';
        slopes.(what).percharac.std    = std(Slopes).';
        slopes.(what).perspan.mean     = mean(transpose(Slopes)).';
        slopes.(what).perspan.std      = std(transpose(Slopes)).';

    elseif strcmp(what,"both")
        slopes.negative = get_slopes(locvals,"negative");
        slopes.positive = get_slopes(locvals,"positive");

        tmp_negative    = slopes.negative.negative;
        slopes          = rmfield(slopes,"negative");
        slopes.negative = tmp_negative;

        tmp_positive    = slopes.positive.positive;
        slopes          = rmfield(slopes,"positive");
        slopes.positive = tmp_positive;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------
function output = replace_zeros(input)

    % if, for any reason, val_min == val_max (too small epsilon for eg)
    % then it becomes and outlier and needs to wether be removed or
    % replaced.
    % Can be by the MEAN or MEDIAN value

    % to find it we take the ceil.
    % as DVALS = VAL_MAXS-VAL_MINS >=0
    % then ceil(DVALS) = 1 if DVALS ~= 0
    % <.>OL stands for OutLier
    if isempty(find(ceil(input)==0,1)) == 0
        if mean(mean(input)) <0
            [xol,yol]   = find(floor(input)==0);
        else
            [xol,yol]   = find(ceil(input)==0);
        end
        % for each column where there are outliers
        for k = 1:length(xol)
            if xol(k)-1 == 0
                input(xol(k),yol(k)) = input(xol(k)+1,yol(k))/2;
            else
                input(xol(k),yol(k)) = mean(input([xol(k)-1,xol(k)+1],yol(k)));
            end
        end
    end

    output = input;