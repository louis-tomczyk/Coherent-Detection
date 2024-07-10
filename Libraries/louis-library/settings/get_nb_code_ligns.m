function count = get_nb_code_ligns(filename)

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

    lines = readlines(filename);
    lines(lines == "" |...
          startsWith(lines,(whitespacePattern|"") + "%"))...
          = [];
    count = numel(lines);


    