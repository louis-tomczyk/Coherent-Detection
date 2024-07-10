function E = Seq2Field(Sequences,varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2024-01-25
%
% ----- Main idea -----%   
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if nargin == 1
        skiprow = 3;
    else
        skiprow = varargin{1};
    end

    shape = size(Sequences);
    if shape(2) > shape(1)
        Sequences = Sequences.';
    end

    Seqs    = Sequences(skiprow:end,:);
    HI      = Sequences(:,2);
    HQ      = Sequences(:,1);
    VI      = Sequences(:,4);
    VQ      = Sequences(:,3);
    
    H       = complex(HI,HQ);
    V       = complex(VI,VQ);
    
    E.field = [H,V];
    E.lambda = 1550;
end