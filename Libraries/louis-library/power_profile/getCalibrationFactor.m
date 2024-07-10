function data = getCalibrationFactor(method)
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET CALIBRATION FACTOR
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-06-20
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Get the calibration factor from exported structures
%
% ----- INPUTS -----
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------

    extractedData   = extractDataFromFiles(method);
    data            = getCalFactor(extractedData,method);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------
% ----- CONTENTS -----
%   extractDataFromFiles
%   getCalFactor
% ---------------------------------------------

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : EXTRACT DATA FROM FILES
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-06-21
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Take data from filename
%
% ----- INPUTS -----
%   DATA_IN     (array) Contains the data obtained with the filenames
%
% ----- BIBLIOGRAPHY -----
%   Functions           :
%   Author              : chat.openai.fr
%   Author contact      : 
%   Date                : 
%   Title of program    : 
%   Code version        : 
%   Type                : 
%   Web Address         : 
% ----------------------------------------------
function extractedData = extractDataFromFiles(method)
    % Sélectionner les fichiers à partir d'une fenêtre de dialogue
    fileNames = uigetfile('*.xml', 'Sélectionner les fichiers', 'MultiSelect', 'on');

    % Si aucun fichier n'est sélectionné
    if isequal(fileNames, 0)
        extractedData = [];
        return;
    end

    % S'assurer que fileNames est une cellule même pour un seul fichier sélectionné
    if ~iscell(fileNames)
        fileNames   = {fileNames};
    end

    % Préallouer le vecteur de résultats
    numFiles        = numel(fileNames);
    extractedData   = zeros(1, numFiles);

    if strcmp(method,'Apeak')
        strBegin = 'Apeak ';
        strEnd = 'calF';
    else
        strBegin = 'calF ';
        strEnd = 'AFC';
    end

    for i = 1:numFiles
        fileName    = fileNames{i};

        % Extraire la donnée numérique entre "Apeak" et "calF" dans le nom du fichier
        startIndex          = strfind(fileName, strBegin) + length(strBegin+' ');
        endIndex            = strfind(fileName, strEnd) - 4;
        dataStr             = fileName(startIndex:endIndex);
        extractedData(i)    = str2double(dataStr);
    end

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET CAL FACTOR
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2023-06-21
%   Version         : 1.1
%
% ----- MAIN IDEA -----
%   Process the peak or calibration factors obtained
%   to return only one value
%
% ----- INPUTS -----
%   DATA_IN     (array) Contains the data obtained with the filenames
%
% ----- BIBLIOGRAPHY -----
% ----------------------------------------------
function data = getCalFactor(data_in,method)

    format long
    
    data        = reshape(data_in,20,[]);
    data_mean   = mean(data);
    data_std    = std(data);
    data        = struct();

    if strcmp(method,'Apeak')
        xdB         = 1:10;
        x           = 10.^(-xdB/10);
        T0          = 1-x;
    
        p           = polyfit(T0,data_mean,1);
        data.cal    = p;
        data.points = [T0.',data_mean.',data_std.'];
    else
        data.points = [data_mean.',data_std.'];
        data.cal    = median(data_mean);
    end