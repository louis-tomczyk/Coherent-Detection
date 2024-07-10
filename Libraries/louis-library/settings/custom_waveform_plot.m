function custom_waveform_plot(fields,axis_struct,n_symbols,file_name_w)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : CUSTOM WAVEFORM PLOT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-30
%   Version         : 1.3
%
% ----- Main idea -----
%   Plot the waveforms at each step of the Digital Signal Processing
%   (DSP) process for polarisation multiplexing or not.
%
% ----- INPUTS -----
%   FIELDS      (structure) with varying size containing all the electric
%                           fields at each DSP step.
%               - IN    (structure) the modulated electric fields
%               - OUT   (structure) the propagated electric fields
%               - CDC   (structure) the CD compensated electric fields
%               - PMDC  (structure) the PMD compensated electric fields
%               - CFOC  (structure) the electric fields with Carrier
%                                   Frequency Offset compensated
%               - CPC   (structure) the electric fields with Carrier
%                                   Phase Compensated
%
%   AXIS_STRUCT (structure) containing at least:
%               - TIME  (array)     Containing the time values
%               - NT    (scalar)    Number of bits per symbol
%
%   N_SYMBOLS   (scalar) The number of samples wanted to be displayed
%   FILE_NAME_W (string) Name of the output figure, string
%
% ----- OUTPUTS -----
%   Figure displayed : A subplot for each DSP step
%   Figure saved with FILE_NAME_W name in current/working directory
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    AA          = 1:n_symbols*axis_struct.Nt;
    XX          = axis_struct.time(AA);

    names_fields= fieldnames(fields);
    n_dsp_steps = length(names_fields);

    if size(fields.in.field,2) == 2*size(fields.in.lambda,1)
        n_polar = 2;
    else
        n_polar = 1;
    end

    f = figure;
    t = tiledlayout(1,n_dsp_steps,'TileSpacing','Compact');
    
    for k=1:n_dsp_steps
        nexttile, % equivalent of subplot, but enables common X/Y labels

            YY.x    = empower(fields.(names_fields{k}).field(:,1));

            if n_polar == 2
                YY.y    = empower(fields.(names_fields{k}).field(:,2));
                YY.xy   = YY.x+YY.y;
                plot(XX,YY.xy(AA,1)','k--')
            end

            hold all
            plot(XX,YY.x(AA,1)','k','linewidth',2)

            if n_polar == 2
                plot(XX,YY.y(AA,1)','color',[0.8,0.8,0.8],'linewidth',2)
            end

            if n_polar == 1
                set_figure_defaults(f,XX,YY,"","","",names_fields{k},9)
            elseif n_polar == 2
                set_figure_defaults(f,XX,YY,"","",["tot","x","y"],names_fields{k},9)
            end
    end
   
    xlabel(t,'time [ns]')
    ylabel(t,'power [mW]')
    exportgraphics(f,file_name_w,'Resolution',100)

end