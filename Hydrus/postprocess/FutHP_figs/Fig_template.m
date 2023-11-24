% Fig_.m
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

matname=fullfile(rootof,'MainFutScenarios.mat');

%% Load pre-proceesed data
load(matname);
disp("Loaded data")
    

%% Make plot



%%
disp("***************************************EOF***************************************")

