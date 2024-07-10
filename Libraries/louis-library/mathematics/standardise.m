function varargout = standardise(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : STANDARDISE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-02-03
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Standardize data by centering by the mean value and normalising the
%   centered input by its standard deviation, with 1/N factor.
%
% ----- INPUTS -----
%   IN      (array)     Signal that we want to standardise
%
% ----- OUTPUTS -----
%   OUT     (array)     Standardized signal
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    varargout = cell(1,3);
    for k = 1:nargin
        in = varargin{k};
        if size(in,1) == 1 && length(size(in)) == 2
            out = (in-mean(in))/std(in,1);
    
        elseif size(in,1) > 1 && length(size(in)) == 2
            out = transpose((in-mean(in))./std(in,1));
            
        elseif length(size(in)) == 3
            out = zeros(size(in,1),size(in,2),size(in,3));
            for k = 1:size(in,3)
                out(:,:,k) = standardise(in(:,:,k));
            end
        end
        varargout{k} = out;
    end

end