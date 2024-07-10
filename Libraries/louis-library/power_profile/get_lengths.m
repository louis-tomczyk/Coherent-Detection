function lengths = get_lengths(locvals,what)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_LENGHTS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-15
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Get the characteristic lengths in the profile from characteristic
%   locations, whether per span or per parameter under study values.
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
    loc_mins            = locvals.loc_mins;
    loc_maxs            = locvals.loc_maxs;
    [Nspans,Ncharac]    = size(loc_mins);

    %%% calculations
    if strcmp(what,"span")
        tmp.maxs            = flipud(loc_maxs);
        tmp.mins            = flipud(loc_mins);

        lengths.(what).maxs = zeros(Nspans-1,Ncharac);
        lengths.(what).mins = zeros(Nspans-1,Ncharac);

        % Npoints gives Npoints-1 intervals
        for k = 1:Nspans-1
            lengths.(what).maxs(k,:)    = tmp.maxs(k,:)-tmp.maxs(k+1,:);
            lengths.(what).mins(k,:)    = tmp.mins(k,:)-tmp.mins(k+1,:);
        end
        
        lengths.(what).percharac.mean   = ((mean(lengths.(what).mins(2:end,:))+mean(lengths.(what).maxs))/2).';
        lengths.(what).percharac.std    = ((std(lengths.(what).mins(2:end,:))+std(lengths.(what).maxs))/2).';

        lengths.(what).perspan.mean     = ((mean(lengths.(what).mins(2:end,:).')+mean(lengths.(what).maxs(2:end,:).'))/2).';
        lengths.(what).perspan.std      = ((std(lengths.(what).mins(2:end,:).')+std(lengths.(what).maxs(2:end,:).'))/2).';

        lengths.(what)  = rmfield(lengths.(what),["mins","maxs"]);

    elseif strcmp(what,"NL")

        lengths.(what).data             = loc_mins(1:end-1,:)-loc_maxs(1:end-1,:);
        lengths.(what).percharac.mean   = mean(lengths.(what).data).';
        lengths.(what).perspan.mean     = mean(lengths.(what).data.').';

        lengths.(what).percharac.std    = std(lengths.(what).data).';
        lengths.(what).perspan.std      = std(lengths.(what).data.').';

        lengths.(what)  = rmfield(lengths.(what),"data");

    elseif strcmp(what,"both")
        lengths.span    = get_lengths(locvals,"span");
        lengths.NL      = get_lengths(locvals,"NL");

        tmp_span        = lengths.span.span;
        lengths         = rmfield(lengths,"span");
        lengths.span    = tmp_span;

        tmp_NL          = lengths.NL.NL;
        lengths         = rmfield(lengths,"NL");
        lengths.NL      = tmp_NL;
    end
end