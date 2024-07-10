function [filenames,Counts] = get_nb_code_ligns_folder(pathh)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_NB_CODE_LIGNS_FOLDER
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-04-02
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get the number of ligns of code in a file
%
% ----- INPUTS -----
%   PATHH       (string)    Folder path where are the codes
%
% ----- OUTPUTS -----
%   COUNTS      (scalar)    Number of ligns of code.
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    cd(pathh)
    nfiles          = length(dir(pathh))-2;
    folder_struct   = dir(pwd);
    Counts          = zeros(1,nfiles);
    filenames       = string(nfiles);

    for k = 1:nfiles
        filename    = folder_struct(k+2).name;
        filenames(k)= filename;
        Counts(k)   = get_nb_code_ligns_int(filename);
    end
    
    filenames   = filenames';
    Counts      = Counts';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   correct_pp
% ---------------------------------------------


% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_NB_CODE_LIGNS
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-04-02
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Get the number of ligns of code in a file
%
% ----- INPUTS -----
%   FILENAME    (string)    Name of the image file of which you want the
%                           Fourier transform.
%
% ----- OUTPUTS -----
%   COUNT       (scalar)    Number of ligns of code.
%
% ----- BIBLIOGRAPHY -----
%   Functions           :
%   Author              : Ive J
%   Author contact      : MathWorks
%   Date                : 2023-04-02
%   Title of program    : 
%   Code version        : 
%   Type                : 
%   Web Address         : https://ch.mathworks.com/matlabcentral/answers/
%                         1939574-how-to-count-uncommented-ligns-in-matlab
%                         -file?s_tid=mlc_ans_email_ques
% ----------------------------------------------

function count = get_nb_code_ligns_int(filename)

    lines = readlines(filename);
    lines(lines == "" |...
          startsWith(lines,(whitespacePattern|"") + "%"))...
          = [];
    count = numel(lines);


    