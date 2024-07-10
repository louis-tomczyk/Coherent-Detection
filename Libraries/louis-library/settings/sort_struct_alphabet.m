function myStruct = sort_struct_alphabet(myStruct)

% Guillaume
% https://fr.mathworks.com/matlabcentral/answers/341950-sort-fieldnames-in-a-structure-alphabetically-ignoring-case
% seen: 2023-02-16 @4pm

    [~, neworder] = sort(lower(fieldnames(myStruct)));
    myStruct = orderfields(myStruct, neworder);