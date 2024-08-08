%
clc;
clear all; 
close all;

list_link_thresh = linspace(0,1,5);
addpath('../videos/')


for index=1:length(list_link_thresh)
    linkThresh = list_link_thresh(index);

    check_only_one_config_no_messenger;

    beep;
    close all;

    figure(101);
    plot_trace_low_graphics();

    beep;

    figure(1);
    res_video_maker_paper();

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
