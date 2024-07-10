function  [locvals,charac] = get_locvals(folder_struct,nfiles,what_charac,charac)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : 
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-21
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Get the locations and values of the extrema in the profile.
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

        % get the number of spans and number of values of the parameter studied.
        % (3) because (1,2) are respectvily './' and '../'.
        [~,shape]           = get_data(folder_struct(3).name);
        Nspans              = shape(1);

        locvals.loc_mins    = zeros(shape);
        locvals.loc_maxs    = zeros(shape);
        locvals.val_mins    = zeros(shape);
        locvals.val_maxs    = zeros(shape);
        
        for k = 1:nfiles
        
            filename    = folder_struct(k+2).name;
        
            if strcmp(what_charac,"epsilon") == 1
                offsets     = [6,4,2];
                tmp_comp    = ["eps_","a"];    
            elseif strcmp(what_charac,"pd") == 1
                offsets     = [3,3,1];          
                tmp_comp    = ["pd_","_"];
            end
            
            index_tmp       = strfind(filename,tmp_comp(1));
            index           = (index_tmp:index_tmp+offsets(1))+offsets(2);
    
            if strcmp(filename(index(end)),tmp_comp(2))
                charac(k)  = str2double(filename(index(1:end-offsets(3))));
            else
                charac(k)  = str2double(filename(index));
            end
    
            %                                      rows   cols
            locvals.loc_mins(:,k)   = get_data(filename,(1:Nspans),1);
            locvals.loc_maxs(:,k)   = get_data(filename,(1:Nspans),2);
            locvals.val_mins(:,k)   = get_data(filename,(1:Nspans),3);
            locvals.val_maxs(:,k)   = get_data(filename,(1:Nspans),4);
        
        end

end