clc; 
close all;


addpath("../functions/")

f = figure(1);
f.Position = [248 396 641 523];

%%% %% Plot single
file_str = "single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat";
% file_str = "single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat";
file_str = "single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_21__i_p2m_30.mat";
load(file_str)
% 
time = size(posArr,5);

tmp = posArr(1,1,:,:,time,1);
pos = reshape(tmp,[2,size(tmp,4)]);

tmp = zpArr(1,1,:,time,1);
z = reshape(tmp,[1,size(tmp,3)]);

tmp = stateArr(1,1,:,time,1);
state = reshape(tmp,[1,size(tmp,3)]);

plot_pos_func(pos, z, state, xx, yy, zz, zEnv, NPop, ArenaScale, linkThresh)


%%% %% Plot All Configs
% load /home/mohsen/Project/Publications/3__messengers/data/large_parameter/MCs_12_Jan_2022_PC_1.0x1.0Arena__50k_tf__lowRand4Explt_cone_BasicMarkov_initENumMsngr_LargeParReg.mat

% N_p2m = length(p2msngrArr);
% N_p2e = length(p2expltArr);
% nSkip = 1;
% 
% n_plts_x = floor(N_p2e/nSkip);
% n_plts_y = floor(N_p2m/nSkip);
% 
% i_MC = 1;
% time = size(posArr,5);
% 
% 

% 
% 
% plt_ctr = 0;
% i_MC = 2;
% for i_p2m=1:nSkip:N_p2m
%     for i_p2e=1:nSkip:N_p2e
%         plt_ctr = plt_ctr + 1;
%         % close all
%         % f = figure(1);
%         % f.Position = [248 396 641 523];
%         hold off
% 
%         % ax = subplot(n_plts_x, n_plts_y, plt_ctr);
% 
%         tmp = posArr(i_p2e,i_p2m,:,:,time,i_MC);
%         pos = reshape(tmp,[2,size(tmp,4)]);
% 
%         tmp = zpArr(i_p2e,i_p2m,:,time,i_MC);
%         z = reshape(tmp,[1,size(tmp,3)]);
% 
%         tmp = stateArr(i_p2e,i_p2m,:,time,i_MC);
%         state = reshape(tmp,[1,size(tmp,3)]);
% 
%         plot_pos_func(pos, z, state, xx, yy, zz, zEnv, NPop, ArenaScale, linkThresh)
% 
%         title_str = sprintf("log_{10} p_e: %1.2f, p_m: %1.2f", log10(p2expltArr(i_p2e)),log10(p2msngrArr(i_p2m)));
% 
%         title(title_str)
% 
%         file_str = sprintf("../codes_for_paper/final_pos_figures/i_p2e__%i__i_p2m__%i",i_p2e,i_p2m);
%         exportgraphics(f,file_str+".png",'resolution',150)
% 
%         % drawnow()
%         % pause(1)
% 
% 
%     end
% end


cbar = colorbar;
caxis([0 1]);
cbar.Limits = [0, 1];
cbar.Label.String = 'Opinion';
set(gca,'FontSize',18)
exportgraphics(f,file_str+"_with_CMap.png",'resolution',150)