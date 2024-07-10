function phi = get_phi(varargin)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_PHI
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-10
%   Version         : 1.3
%
% ----- Main idea -----
%   Get the Phase of the input field
%
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to be normlised
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
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

    % if NARGIN == 1 then no need to downsample
    % if NARGIN == 2 then need for downsampling at 1SPS

    Einput  = varargin{1};
    n_polar = size(Einput.field,2);

    if nargin == 1
        assert(isstruct(varargin{1}) == 1,"Argument should be EINPUT")
        for k = 1:n_polar
            phi.(strcat("pol_",num2str(k))) = angle(Einput.field(:,k));
        end
    elseif nargin == 2
        assert(length(fieldnames(varargin{1})) == 2,"Argument should be EINPUT")
        assert(length(fieldnames(varargin{2})) > 15,"Argument should be AXIS")

        Axis = varargin{2};
        if size(Einput.field,1) > Axis.Nsymb
            Einds = ds(Einput,1,Axis.Nsps,1);
            for k = 1:n_polar
                phi.(strcat("pol_",num2str(k))) = angle(Einds.field(:,k));
            end
        else
            for k = 1:n_polar
                phi.(strcat("pol_",num2str(k))) = angle(Einput.field(:,k));
            end
        end

    end

    
