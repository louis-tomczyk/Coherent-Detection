function init_step(varargin)

    % ---------------------------------------------
    % ----- INFORMATIONS -----
    %   Function name   : INIT_STEP
    %   Author          : louis tomczyk
    %   Institution     : Telecom Paris
    %   Email           : louis.tomczyk@telecom-paris.fr
    %   Date            : 2023-11-03
    %   Version         : 1.2.7
    %
    % ----- Main idea -----
    %   Initialises the paths access to functions and clears all
    %
    % ----- INPUTS -----
    %   WHAT    (string)    Can whether be "nli", "PPE" or nothing
    %
    % ----- OUTPUTS -----
    % ----- BIBLIOGRAPHY -----
    % ---------------------------------------------
    
        clearvars -except varargin
        close all
        clc
           

        if nargin == 0
            cd(pwd)  
        elseif nargin > 1
            if strcmp(varargin{1},'local')
                cd('/home/louis/Documents/6_Telecom_Paris/3_Codes/0_louis/1_Optilux/')
            end
    
            what = varargin{2};
            if isStringInStringArray(what,"nli") == 1
                addpath(strcat(pwd,"/Libraries/louis library/anli"))
            elseif isStringInStringArray(what,"PPE") == 1
                addpath(pwd+"/Libraries/louis library/power_profile/")
            end
        else
            what = varargin{1};
            cd('/home/louis/Documents/6_Telecom_Paris/3_Codes/0_louis/1_Optilux/')  
            if isStringInStringArray(what,"nli") == 1
                addpath(strcat(pwd,"/Libraries/louis library/anli"))
            elseif isStringInStringArray(what,"PPE") == 1
                addpath(pwd+"/Libraries/louis library/power_profile/")
            end
        end
        
        addpath("./")
        addpath("./Libraries/Elie library")
        addpath("./Libraries/Optilux library")
        addpath("./Libraries/louis library")
        addpath("./Libraries/louis library/classic_dsp/")
        addpath("./Libraries/louis library/constants/")
        addpath("./Libraries/louis library/mathematics/")
        addpath("./Libraries/louis library/propagation/")
        addpath("./Libraries/louis library/settings/")
        addpath("./Libraries/louis library/transmitter/")
    
        format long
    
        global Axis
    
        set(0,'defaultfigurecolor',[1 1 1])
        set(groot,'defaultAxesTickLabelInterpreter','latex'); 
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        set(groot, 'defaultFigureUnits','normalized')
        set(groot, 'defaultFigurePosition',[0 0 1 1])


function bool = isStringInStringArray(stringArray,myString)
    % ---------------------------------------------
    % ----- INFORMATIONS -----
    %   Function name      IS STRING IN STRINGARRAY
    %   Author             louis tomczyk
    %   Institution        Telecom Paris
    %   Email              louis.tomczyk@telecom-paris.fr
    %   Date               2023-02-24
    %   Version            1.0
    %
    % ----- Main idea -----
    % ----- INPUTS -----
    % ----- OUTPUTS -----
    % ----- BIBLIOGRAPHY -----
    % ---------------------------------------------
    
        assert(nargin == 2,"too much or too few arguments")
        assert(class(stringArray)=="char" || class(stringArray)=="string",...
            "arg 1 should be a string")
        
        bool = isempty(find(ismember(stringArray,myString), 1)) == 0;