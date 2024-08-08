%
clc;
clear all; close all;

folder_path = 'data/';
files = dir(fullfile(folder_path, '*.mat'));


parfor k = 1:length(files)

% for k = 3:3 % length(files):length(files) % 3:4
    filename = fullfile(folder_path, files(k).name);
    fprintf('Processing %s\n', filename);
    
    % Load the .mat file data (optional, depends on needs)
    % data = load(filename);

    opinion_plot_func(filename);

    beep;
    pause(2);
    beep;
end


beep;
beep;
pause(2);
beep;
pause(1);
beep
