%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /home/alumno.upv.es/almousa/CMCP/Trabajo/plotsmatlab/Medidas Trabajo CMCP.xlsx
%    Worksheet: Plots
%
% Auto-generated by MATLAB on 17-Nov-2021 17:26:36

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = "Plots";
opts.DataRange = "B5:O12";

% Specify column names and types
opts.VariableNames = ["Hilos", "TiempoFuerte400x400", "TiempoFuerte1000x1000", "TiempoDebil", "SpeedupFuerte400x400", "SpeedupFuerte1000x1000", "SpeedupDebil", "CosteFuerte400x400", "OverheadFuerte400x400", "CosteFuerte1000x1000", "OverheadFuerte1000x1000", "EficienciaFuerte400x400", "EficienciaFuerte1000x1000","EficienciaDebil"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
OpenMP = readtable("/home/alumno.upv.es/almousa/CMCP/Trabajo/plotsmatlab/Medidas Trabajo CMCP.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = "Plots";
opts.DataRange = "B17:O24";

% Specify column names and types
opts.VariableNames = ["Procesos", "TiempoFuerte400x400", "TiempoFuerte1000x1000", "TiempoDebil", "SpeedupFuerte400x400", "SpeedupFuerte1000x1000", "SpeedupDebil", "CosteFuerte400x400", "OverheadFuerte400x400", "CosteFuerte1000x1000", "OverheadFuerte1000x1000", "EficienciaFuerte400x400", "EficienciaFuerte1000x1000","EficienciaDebil"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
MPI = readtable("/home/alumno.upv.es/almousa/CMCP/Trabajo/plotsmatlab/Medidas Trabajo CMCP.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

