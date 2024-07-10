function distance = get_distance(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_DISTANCE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-17
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   We aim to get the distance between an INPUT and a REFERENCE object in a
%   N dimensions space with a given (or not) calculation METHOD.
%   Default method is L² distance.
%
% ----- INPUTS -----
%   varargin{1}-INPUT       (array)
%   varargin{2}-REFERENCE   (array)
%   varargin{3}-METHOD{optional}
%                           (string)    Choose among L1 - L2
% ----- OUTPUTS -----
%   DISTANCE    (array) distances between the objects and the references
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------


    input       = varargin{1};
    reference   = varargin{2};

    if nargin == 3
        method  = varargin{3};
    else
        method  = 'L2';
    end 

    [nr,nc]     = size(input);
    distance    = zeros(nr,nc);

    for k = 1:nr
        for j = 1:nc
            if strcmp(method,"L1") == 1
                distance(k,j) = abs(input(k,j)-reference);
            elseif strcmp(method,'L2') == 1
                distance(k,j) = (abs(input(k,j)-reference)).^2; 
            end
        end
    end

end