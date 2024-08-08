%
clc;
clear all; close all;

addpath("../functions/")

% list_i_p2e = [01, 49, 07, 21, 36, 15, 06];
% list_i_p2m = [49, 01, 12, 30, 45, 02, 07];

% list_i_p2e = [07, 05, 11, 12, 11, 12, 07, 08, 07, 35, 36, 36];
% list_i_p2m = [08, 06, 08, 08, 07, 07, 14, 15, 13, 45, 46, 47];

% list_i_p2e = [44:-1:39, 44:-1:39, 44:-1:39];
% list_i_p2m = [49*ones(1,6), 48*ones(1,6), 47*ones(1,6)]; % 49*ones(size(list_i_p2e)); 
% 
% 
% parfor index=1:length(list_i_p2e)
%     i_p2e = list_i_p2e(index);
%     i_p2m = list_i_p2m(index);
% 
%     % check_only_one_config;
%     strng = check_only_one_config_func(i_p2e, i_p2m);
% 
%     beep;
% 
%     % figure(101);
%     % plot_trace_low_graphics();
% 
%     beep;
% 
%     figure(1);
%     % res_video_maker_paper();
%     % res_video_maker_paper_with_opinion_distr();clos
%     res_video_maker_paper_with_opinion_distr_and_plot(strng);
% 
%     beep;
%     pause(2);
%     beep;
% end

folder_path = 'selected_data/';
files = dir(fullfile(folder_path, '*.mat'));

% parfor k = 1:length(files)

% for k = 3:3 % length(files):length(files) % 3:4
    % filename = fullfile(folder_path, files(k).name);
    % fprintf('Processing %s\n', filename);

    filename = "selected_data/single_run_31_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_40__i_p2m_49.mat";
    
    % Load the .mat file data (optional, depends on needs)
    % data = load(filename);

    res_video_maker_paper_with_opinion_distr_and_plot(filename);

    beep;
    pause(2);
    beep;
% end

beep;
pause(2);
beep;
pause(1);
beep
