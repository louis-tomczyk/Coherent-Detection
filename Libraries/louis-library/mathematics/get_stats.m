function out = get_stats(in,what,q)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET STATS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-08
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get statistics from signal
%
% ----- INPUTS -----
%   IN      (structure/array)   Signal from which we want the statistics
%   WHAT    (string)            What (empirical) statistics we want, 
%                               choose among: 
%                               MEAN    - mean
%                               VAR     - VARiance with 1/N factor
%                               STD     - STandard Deviation with 1/N
%                                       factor
%                               SKW     - SKeWness
%                               KTS     - KurToSis
% ----- OUTPUTS -----
%   OUT      (structure/array)  The statistics we wanted.
%                               Array case if: 
%                               - real input - 
%                               - complex input but with real or imaginary 
%                                   parts exclusively
%                               Structure if both real/imaginary parts
%                               wanted
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    if isstruct(in) == 1
        in = in.field;
    end
    
    if isreal(in)
        switch what
            case "mean"
                out = mean(in);
            case "var"
                out = var(in,1);
            case "std"
                out = std(in,1);
            case "skw"
                out = skewness(in);
            case "kts"
                out = kurtosis(in);
            otherwise
                what= input("not known, please select among 'mean','var','std','skw', 'kts'     ","s");
                out = get_stats(in,what);
        end
    else
        if q == 'r'
            switch string(what)
                case "mean"
                    out= mean(real(in));
                case "var"
                    out= var(real(in),1);
                case "std"
                    out= std(real(in),1);
                case "skw"
                    out= skewness(real(in));
                case "kts"
                    out= kurtosis(real(in));
                otherwise
                    what    = input("not known, please select among 'mean','var','std','skw', 'kts'     ","s");
                    out     = get_stats(in,what);
            end
        elseif q == 'i'
            switch string(what)
                case "mean"
                    out = mean(imag(in));
                case "var"
                    out = var(imag(in),1);
                case "std"
                    out = std(imag(in),1);                
                case "skw"
                    out = skewness(imag(in));
                case "kts"
                    out = kurtosis(imag(in));                
                otherwise
                    what    = input("not known, please select among 'mean','var','std','skw', 'kts'     ","s");
                    out     = get_stats(in,what);
            end
        elseif q == 'b'
            switch string(what)
                case "mean"
                    out.r= mean(real(in));
                    out.i = mean(imag(in));
                case "var"
                    out.r= var(real(in),1);
                    out.i = var(imag(in),1);
                case "std"
                    out.r= std(real(in),1);
                    out.i = std(imag(in),1);                
                case "skw"
                    out.r= skewness(real(in));
                    out.i = skewness(imag(in));
                case "kts"
                    out.r= kurtosis(real(in));
                    out.i = kurtosis(imag(in));                
                otherwise
                    what    = input("not known, please select among 'mean','var','std','skw', 'kts'     ","s");
                    out     = get_stats(in,what);
            end
        end
    end

end