% addpath( genpath( '../code/' ) );
% addpath( genpath( '../data/' ) );
% addpath( genpath( '../external/' ) );
% set(0,'DefaultFigureWindowStyle','docked')
restoredefaultpath
clear RESTOREDEFAULTPATH_EXECUTED

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

clear all; close all; clc;
