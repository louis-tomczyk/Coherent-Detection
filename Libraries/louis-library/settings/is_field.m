function isFieldResult = is_field (inStruct, fieldName)
   

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : IS_FIELD
%   Author          : MathWorks SUpport team
%   Institution     : MathWokds
%   Email           : 
%   Date            : 2017-06-12
%   Version         : 1.0
%
% ----- MAIN IDEA -----
%   Find if field in ANYwhere in a structure (nested levels included)
% ----- INPUTS -----
%   INSTRUCT    (string/array)  name of the structure or an array of
%                               structures to search
%   FIELDNAME   (string)        name of the field for which the function
%                               searches
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    
    isFieldResult = 0;
    f = fieldnames(inStruct(1));
    for i=1:length(f)
        if(strcmp(f{i},strtrim(fieldName)))
            isFieldResult = 1;
            return;
        elseif isstruct(inStruct(1).(f{i}))
            isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
        end

        if isFieldResult
            return;
        end
    end
end