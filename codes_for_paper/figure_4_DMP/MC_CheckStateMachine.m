% Title: Run A MonteCarlo study to evaluate Dichotomous Markov Process
% Description: This script runs MoteCarlo Simulations
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: May 27, 2024
% Version: 2.0
% Usage: Run this script with MATLAB R2023a or later.
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"

clc;
clear all;
close all;

tf = 50e3;                      % Duration of Experiment
NPop = 100;                     % Population Size

nTVars = 10;                   % Number of Time Samples
nSkipSave = round(tf/nTVars);   % Skip every x steps

% parpool(20)

nMC = 200;                        % Number of MonteCarlo Simulations

p2expltArr = exp(-[1:2:50]./5);     % Array for P_
p2msngrArr = p2expltArr(1:end-1);
np2E = length(p2expltArr);
np2M = length(p2msngrArr);


% expltDur = nan(nMC,NPop,np2E,np2M);
numExpl = nan(nMC,nTVars,np2E,np2M);
% numEdges = nan(nMC,NPop,np2E,np2M);

avgNumExpl = nan(nMC,np2E,np2M);

% 0 : Messenger
% 1 : Exploiter

wbar = waitbar(0, 'Starting');

for ip2E=1:np2E
    p2E = p2expltArr(ip2E);
    waitbar(ip2E/np2E, wbar, sprintf('Progress: %d %%', floor(ip2E/np2E*100)));

    disp(['#P2E: ',num2str(ip2E),' (out of: ',num2str(np2E),'), P2E: ',num2str(p2E)]);

    for ip2M=1:np2M
        p2M = p2msngrArr(ip2M);
        disp(['----- | #P2M: ',num2str(ip2M),' (out of: ',num2str(np2M),'), P2M: ',num2str(p2M)]);

        parfor iMC=1:nMC

%             state = zeros(NPop,1); % start with all Messengers
%             state = ones(NPop,1); % start with all Exploiters

            ExpNExploiter = p2E/(p2M+p2E);
            tmp = rand(NPop,1);
            state = tmp>ExpNExploiter;

            stArr = nan(NPop,tf);

            for t=1:tf
                stArr(:,t) = state;

                randArr = rand(NPop,1);

                posib2msngr = p2M*state;
                posib2explt = p2E*not(state);

                expChange = randArr<posib2msngr;
                msgChange = randArr<posib2explt;

                state(expChange) = not(state(expChange));
                state(msgChange) = not(state(msgChange));
            end

            % expltDur(iMC,:,ip2E,ip2M) = mean(stArr,2);          % ratio of time each agent stayed in Exploiter state

            tmp = mean(stArr);
            numExpl(iMC,:,ip2E,ip2M) = tmp(1:nSkipSave:end);    % mean ratio of pop. doing exploitation

            % tmp = diff(stArr,1,2);
            % numEdges(iMC,:,ip2E,ip2M) = sum(abs(tmp),2);        % number of switches

            avgNumExpl(iMC, ip2E, ip2M) = mean(stArr,'all');

        end
    end
end
% close(wbar)
% to check if it is updated!

% save Res_CheckStateMachine_expInit_2024_highMC_high_T.mat
save Res_CheckStateMachine_expInit_2024_new.mat

% plot_DMP_sojourn_time();
plot_DMP_messengers_ratio();
% cd ../videos/
% run check_and_save_single_runs_bunch.m
