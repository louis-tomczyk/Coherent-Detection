function remove_rows_or_cols(filename,row_or_col,step)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : 
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2022-11-15
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Remove rows or columns of a matrix saved in a file
%
% ----- INPUTS -----
%   FILENAME    (string)    Name of the output file
%   ROW_OR_COL  (string)    What we want to delete
%   STEP        (scalar)    
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

% values for ROW_OR_COL
%   - ROW --- keep one row every 2 rows
%   - COL --- keep one col every 2 cols

% values for ODD_OR_EVEN
%   - ODD --- start from 1st ROW/COL
%   - EVEN --- start from 2nd ROW/COL

data = get_data(filename);

if strcmp(row_or_col,'row') == 1
    data_out = data(2:step:end,:);
elseif strcmp(row_or_col,'col') == 1
    data_out = data(1,2:step:end,:);
end

tmp         = char(filename);
new_name    = strcat(tmp(1:end-4),'_cleaned.txt');
writematrix(data_out,new_name)
