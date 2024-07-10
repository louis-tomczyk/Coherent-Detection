function varargout = qpsk_demod_choice(method,varargin)
    % ---------------------------------------------
    % ----- INFORMATIONS -----
    %   Function name   : QPSK_DEMOD - Quadrature Phase Shift
    %                       Keying Demodulation
    %   Author          : louis tomczyk
    %   Institution     : Telecom Paris
    %   Email           : louis.tomczyk.work@gmail.com
    %   Date            : 2022-09-08
    %   Version         : 1.0
    %
    % ----- Main idea -----
    %   Demodulation of QPSK constellation with data-aided
    %   approach
    % 
    % ----- INPUTS -----
    % ----- OUTPUTS -----
    % ----- BIBLIOGRAPHY -----
    %   Articles
    %   Author              : Nicolas ROSSI, Petros RAMANTANIS, Jean Claude ANTONA
    %   Title               : Nonlinear Interferenc ce Noise Statistics in Unmanagged Coherent
    %                         Networks with Channnels Propagating over Differennt Lightpaths
    %   Jounal              : ECOC
    %   Volume - N°         : NA
    %   Date                : 2014-11-24
    %   DOI                 : 10.1109/ECOC.2014.6964043
    % ---------------------------------------------
    
    %%% MAINTENANCE

    if strcmp(method,"data aided") == 1
        % checking the number of input arguments
        assert(nargin == 2,"Wrong number of input arguments." + ...
            "Should be either (EINPUT,AXIS) or (EINPUT,PHI).")

        % getting the name of input arguments
        argnames = string(nargin);
        for k=1:nargin
            argnames(k)  = inputname(k);
        end
        
        % checking if first argument is the input field
        assert(length(fieldnames(varargin{1})) == 2,"First argument should be EINPUT")
        Einput  = varargin{1};

        % checking what input arguments there are
        % if AXIS is there
        if isempty(find(strcmp(argnames,"Axis"),1)) == 0
            Axis    = varargin{2};
        % if PHI is there
        else
            Phi     = varargin{2};
        end

        n_polar = size(Einput.field,2);

        %%% PHASE VALUES
        % if PHI does not exists, it means that the EINPUT is the modulated
        % signal and then we need the PHI values
        % it PHI exists, it means that the EINPUT is the propagated (or other)
        % signal.
        if exist("Phi",'var') == 0
            Phi     = get_phi(Einput,Axis);
            if size(Einput.field,1) > Axis.Nsymb
                Einput  = ds(Einput,1,Axis.Nsps,1);
            end
        end
        
        %%% DEMODULATION
        Edemod  = Einput;
        for k = 1:n_polar
            Edemod.field(:,k) = (Einput.field(:,k)).*...
                                fastexp(-Phi.(strcat("pol_",num2str(k))));
        end

        %%% OUTPUTS
        if isempty(find(strcmp(argnames,"Axis"),1)) == 0
            varargout{1} = Edemod;
            varargout{2} = Phi;
        else
            varargout{1} = Edemod;
        end
    
    elseif strcmp(method,"distance") == 1
        n_samp  = size(in.field,1);
        n_polar = size(in.field,2);
        out     = in;

        %%% downsampling if necessary
        if n_samp > Axis.Nsymb
            % getting the number of samples per symbol
            Nsps    = n_samp/Axis.Nsymb;
            assert(is_integer(Nsps) == 1,"Number of samples per symbol should be an integer")

            % downsample to get 1 sample per symbol
            in          = ds(in,1,Nsps,1);
            n_samp      = Axis.Nsymb;
            out         = rmfield(out,"field");
            out.field   = zeros(n_samp,n_polar);
        end
        
        m   = [-3,-1,+1,+3];
        phi = exp(1i*m*pi/4);

        %%% demodulation
        for j = 1:n_polar
            for k = 1:n_samp
                amp             = abs(in.field(k,j));
                th              = fastexp(angle(in.field(k,j)));
                distance        = get_distance(phi,th,'L1');
                [~,imin]        = min(distance);
                m_hat           = m(imin);
                out.field(k,j)  = amp.*th*fastexp(-m_hat*pi/4);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------
function phi = get_phi(varargin)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_ANLI
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

    

