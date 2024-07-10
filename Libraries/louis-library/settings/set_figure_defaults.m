function set_figure_defaults(f,XX,YY,xlab, ylab, leg, t,fsize )

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_FIGURE_DEFAULT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-23
%   Version         : 1.0
%
% ----- Main idea -----
%   Set default figure parameters.
%   
% ----- INPUTS -----
%   F   (matlab.ui)     the figure
%   XX  (vector)        the x-axis inputs
%   YY  (structure)     the y-axis inputs
%           - YY.x the x-polarisation field
%           - YY.y the y-polarisation field (optional)
%           - YY.x the x and y-polarisation field (optional but mandatory
%           if YY.y exists)
%   XLAB    (string)    label for x-axis
%   YLAB    (string)    label for y-axis
%   LEG     (array)     legend entries
%   T       (string)    title
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    % DOMAINS

    N = size(fieldnames(YY),1);
    assert(N==1|N==3,"     Signal must have 1 or 2 polarisation");

    if N == 1
        % YY.x
        n_polar = 1;
    else
        % YY.x
        % YY.y
        % YY.xy
        n_polar = 2;
    end

    scale = 1.25;
    if n_polar == 2
        axis([0.75*min(XX),1.1*max(XX),...
            0.25*min(min([YY.x,YY.y,YY.xy])),scale*max(max([YY.x,YY.y,YY.xy]))]);
    else
        axis([0.75*min(XX),scale*max(XX),...
            0.25*min(YY.x),scale*max(YY.x)]);   
    end

    assert(fsize-5>2,"\n\t Fontsize must be >7");
    % GCA
    set(gca,"fontsize",fsize-4);
    
    % XLABS
    xlabel(xlab,"fontsize",fsize-2,"fontweight","bold");

    % YLABS
    ylabel(ylab,"fontsize",fsize-2,"fontweight","bold");

    % LEGEND
    if strcmp(leg,"") == 0
        if n_polar == 2
            legend(leg,"fontsize",fsize-2,"fontweight","bold",'Orientation','vertical','color','none');
            legend boxoff
        end
    end
    % TITLE
    title(t,"fontsize",fsize,"fontweight","bold");
end