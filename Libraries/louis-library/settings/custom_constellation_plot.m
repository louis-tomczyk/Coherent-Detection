function custom_constellation_plot(fields,file_name_w)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CUSTOM CONSTELLATION PLOT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-07-13
%   Version         : 1.4
%
% ----- Main idea -----
%   Plot the constellations at each step of the Digital Signal Processing
%   (DSP) process for polarisation multiplexing or not.
%
% ----- INPUTS -----
% if Coherent Detection only
%   FIELDS  (structure) with varying size containing all the electric fields 
%           at each DSP step.
%               - IN    : the modulated electric fields
%               - OUT   : the propagated electric fields (optional)
%               - CDC   : the CD compensated electric fields (optional)
%               - PMDC  : the PMD compensated electric fields (optional)
%               - CFOC  : the electric fields with Carrier Frequency Offset
%                           compensated (optional)
%               - CPC   : the electric fields with Carrier Phase
%                           Compensated (optional)
% 
% else if Power Profile Estimation
%   FIELDS  (structure) with varying size containing all the electric fields 
%           at each DSP step.
%               - IN        : the modulated electric fields
%               - PD        : the pre-dispersed field before channel
%               - OUTLINK   : the field after propagation into the whole
%                               link
%               - CDC       : the CD compensated electric fields
%               - INPPE     : the electric fields before the PPE
%
%   FILE_NAME_W: name of the output figure, string
%
% ----- OUTPUTS -----
%   Figure displayed : 1st row --- X-pol constellations
%                      2nd row --- Y-pol constellations
%   Figure saved with FILE_NAME_C name in current/working directory
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

%% maintenance
    % new file name creation
    file_name_c = strrep(file_name_w,"waveform","constellation");

    % get the names of all the input fields corresponding to the different
    % DSP steps
    names       = fieldnames(fields);

    % number of DSP steps and polarisation to determine the number of rows
    % and columns for subplot.

    if isfield(fields,'norm') == 1
        % '-1' as the first field is the boolean to normalise or not
        n_dsp_steps = size(names,1)-1;
    else
        fields.norm = 0;
        n_dsp_steps = size(names,1);
    end
    n_polars    = size(fields.(names{1}).field,2);

    % we want to normalise the constellation to put it in a SAME square of
    % size length equals to 2 for ALL the DSP steps for ALL polarisation states
    %   -1 <= REAL(fields) <= 1
    %   -1 <= IMAG(fields) <= 1

    % initialisation of the maxima
%     max_R   = -Inf*ones(1,n_polars);
%     max_I   = -Inf*ones(1,n_polars);

    max_R   = zeros(n_polars,n_dsp_steps);
    max_I   = zeros(n_polars,n_dsp_steps);
    rho     = zeros(n_polars,n_dsp_steps);

    % update the maxima for each DSP step
    for k=1:n_dsp_steps

        % create a temporary variable to not touch the input structure
        if fields.norm == 1
            names_norm  = strcat(names{k},'_norm');
        end
% for each polarisation state we compare the previous maximum with
% the new maximum and take the highest value
%         for j=1:n_polars
%             tmp     = fields.(names{k}).field(:,j);
%             max_R(j)= max([max_R(j), max(real(tmp))]);
%             max_I(j)= max([max_I(j), max(imag(tmp))]);
%         end
% 
%         % radius calculation for normalisation by it
%         rho = abs(max(max_R)+1i*max(max_I));
%         fields.(names_norm).field = fields.(names{k}).field_ds/rho;

%% scaling
% we scale the constellation according to each polarisation state and DSP
% step
        if fields.norm == 1
            for j=1:n_polars
                tmp = fields.(names{k}).field(:,j);
                max_R(j,k)= max(real(tmp));
                max_I(j,k)= max(imag(tmp));
    
                % radius calculation for normalisation by it
                rho(j,k) = abs(max_R(j,k)+1i*max_I(j,k));
                fields.(names_norm).field(:,j) = fields.(names{k}).field_ds(:,j)/rho(j,k);
            end
        end

    end

%% plotting
    f = figure;
    t = tiledlayout(n_polars,n_dsp_steps);

    % plotting for each polarisation
    for k=1:n_dsp_steps
        nexttile,
        
            if isfield(fields,'norm') && fields.norm == 1
                names_norm  = strcat(names{k},'_norm');
                scatter(real(fields.(names_norm).field(:,1)),...
                    imag(fields.(names_norm).field(:,1)),...
                    'MarkerEdgeColor','k','MarkerFaceColor','k')
    
                % for clarity purpose enlarge a bit the square limits 
                % just in case
                axis([-1.1,1.1,-1.1,1.1])
            else
                scatter(real(fields.(names{k}).field(:,1)),...
                    imag(fields.(names{k}).field(:,1)),...
                    'MarkerEdgeColor','k','MarkerFaceColor','k')

            end

%             set(gca,'xticklabel',[])
%             set(gca,'yticklabel',[])
            title(sprintf("%s - x",names{k}))

            % square ratio of each subplots
            pbaspect([1,1,1])
        if n_polars == 2
            nexttile(k+n_dsp_steps),
                scatter(real(fields.(names_norm).field(:,2)),...
                        imag(fields.(names_norm).field(:,2)),...
                        'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8])
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                title(sprintf("%s - y",names{k}))
                axis([-1.1,1.1,-1.1,1.1])
                pbaspect([1,1,1])
        end
    end
    xlabel(t,'In Phase [normalised]')
    ylabel(t,'Quadrature [normalised]')

    % saving the figure
    exportgraphics(f,file_name_c,'Resolution',600)

end