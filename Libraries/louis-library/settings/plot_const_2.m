function t = plot_const_2(filename,varargin)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PLOT_CONST
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-08-13
%   Version         : 1.3
%
% ----- Main idea -----
%   Plot the constellations
%
% ----- INPUTS -----
%   VARARGIN{k}:(structure) containing the Field - MANDATORY
%               - LAMBDA [nm]: wavelength
%               - FIELD [sqrt(mW)]: normalised electric fields
%  for k = 1:nargin
%
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    fig = figure;
    tmp = 0;
    for k = 1:nargin-1
        tmp = tmp + size(varargin{k}.field,2);
    end

    if tmp == 1
        t = tiledlayout(1,1);
    elseif tmp == 2
        t = tiledlayout(1,tmp);
    
    else
        t = tiledlayout(2,ceil(tmp/2),'TileSpacing','Compact');
    end
    
    for k = 1:nargin-1
        nexttile,

        if size(varargin{k},2) == 1
            Etmp = varargin{k};
            if isstruct(Etmp) == 1
                if isfield(Etmp,'field_ds') == 1
                    E = Etmp.field_ds;
                else
                    E = Etmp.field;
                end
            else
                E = Etmp;
            end
            
            if size(E,2) == 1
                I = real(E);
                Q = imag(E);
    
                scatter(I,Q,20,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')

                xx  = 1.1*max(abs(I));
                yy  = 1.1*max(abs(Q));
                z   = max(max(xx,yy));
                axis([-1,1,-1,1]*z)

                %title(inputname(k),'Color','k')
            else
                Ex = E(:,1);
                Ey = E(:,2);

                Ix = real(Ex);
                Iy = real(Ey);

                Qx = imag(Ex);
                Qy = imag(Ey);
    
                scatter(Ix,Qx,20,'filled','MarkerEdgeColor','y','MarkerFaceColor','y')

                xx  = 1.1*max(abs(Ix));
                yy  = 1.1*max(abs(Qx));
                z   = max(xx,yy);
                axis([-1,1,-1,1]*z)

                title(strcat(inputname(k),' X'),'Color','w')
                pbaspect([1,1,1])
                set(gca,'color','k',...
                    'Xcolor','w','Ycolor','w')
                set(gcf,'color','k')

                %%%%%%%%%%%%%%%%%%%%%%%%%
                nexttile,
                scatter(Iy,Qy,20,'filled','MarkerEdgeColor','y','MarkerFaceColor','y')

                xx  = 1.1*max(abs(Iy));
                yy  = 1.1*max(abs(Qy));
                z   = max(xx,yy);
                axis([-1,1,-1,1]*z)

                %title(strcat(inputname(k),' Y'),'Color','w')

                set(gca,'color','k',...
                    'Xcolor','w','Ycolor','w')
                set(gcf,'color','k')
            end
        else
            colors  = ['y','b','r','m'];
            Etmps   = varargin{k};
            N       = size(Etmps,2);

            assert(N<=length(colors),'too many constellations to surimpose')
            hold on

            xx = 0;
            yy = 0;
            for j = 1:N
                Etmp = Etmps(j);

                if isstruct(Etmp) == 1
                    E = Etmp.field_ds;
                else
                    E = Etmp;
                end
                
                I = real(E);
                Q = imag(E);
                scatter(I,Q,20,'filled','MarkerEdgeColor',colors(j),'MarkerFaceColor',colors(j))
                xx  = 1.1*max(xx,max(abs(I)));
                yy  = 1.1*max(yy,max(abs(Q)));
                
            end

            z = max(xx,yy);
            axis([-1,1,-1,1]*z)

        end

        pbaspect([1,1,1])
    
        % change w for k and k for w
        set(gca,'color','w',...
            'Xcolor','w','Ycolor','w')  
        set(gcf,'color','w')
    end

    if tmp > 1
        xlabel(t,'I','color','w')
        ylabel(t,'Q','color','w')
    else
        xlabel('I','color','w')
        ylabel('Q','color','w')
    end

    exportgraphics(fig,strcat(filename,'.png'),"Resolution",300)
end