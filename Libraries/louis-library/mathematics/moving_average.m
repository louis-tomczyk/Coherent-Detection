function vector_out = moving_average(vector_in,params)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : MOVING_AVERAGE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-05-23
%   Version         : 1.0
%
% ----- Main idea -----
%   Average a vector by moving blocks
%
% ----- INPUTS -----
%   VECTOR_IN   (aray)      the vector to average
%   PARAMS      (structure) the averaging parameters structure
%                   - PERIOD    How many points used for averaging
%                   - METHOD    Two methods are implented
%                       - "mirror" : consists in repeating the vector by
%                       concatenating the flipped version of the initial
%                       vector. For examples:
%                           vector_in = [1,2,3] 
%                           ---> vector_in_bis = [1,2,3,3,2,1]
%                       Both output and input have the same lengths.
%     
%                       - "no fake": consists in not repeating the vector
%                       which does not add "fake values".
%                       length(vector_out) = length(vector_in)-PERIOD
%                  
% ----- OUTPUTS -----
%  VECTOR_OUT : The average vector which size depends on the method.
%               For details, look at PARAMS.METHOD
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    % number of elements in the vector to average
    N = length(vector_in);

    % number of samples used for averaging
    M = params.av_period;

    if M == 1
        vector_out = vector_in;
        return
    end

    if strcmp(params.av_method,'mirror')

        % make a copy of the original vector, flip from LEFT to RIGHT the
        % copy and concatenate along the 2nd axis (columns) both
        vector_in_bis = [vector_in,fliplr(vector_in)];

        % pre-allocate memory space
        vector_out = zeros(1,N);

        % average the original vector over the period and add its value to
        % the new vector
        for k=1:length(vector_in)
            vector_out(k) = 1/M*sum(vector_in_bis(k:k+M-1));
        end
        
    % if other string is used for the method, implentation of bloc
    % averaging without adding data.
    else
        vector_out = zeros(1,N-M);
        for k=1:length(vector_in)-M
            vector_out(k) = 1/M*sum(vector_in(k:k+M-1));
        end
    end
    
%     disp("averaging finished")