% Title: Running a single simulation of the DMP (state machine) for a
% specific parameter of the DMP
% Description: This script runs a sample of DMP and visualizes data for
% checkup
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: May 27, 2024
% Version: 2.0
% Usage: Run this script with MATLAB R2023a or later. 
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code 
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"

clc; clear all; close all;
tf = 10e3;
NPop = 100;100;

nTVars = 200;
nSkipSave = round(tf/nTVars);


nMC = 1;

% stMat = nan(NPop,nTVars,nMC);


p2expltArr = exp(-[1:2:100]./5);
p2msngrArr = p2expltArr(1:end-1);
np2E = length(p2expltArr);
np2M = length(p2msngrArr);


expltDur = nan(nMC,NPop,np2E,np2M);
numExpl = nan(nMC,nTVars,np2E,np2M);
numEdges = nan(nMC,NPop,np2E,np2M);

for ip20=15%1:np2E
    p20 = p2expltArr(ip20);
    disp(['#P2E: ',num2str(ip20),' (out of: ',num2str(np2E),'), P2E: ',num2str(p20)]);
    
    for ip21=20%1:np2M
        p21 = p2msngrArr(ip21);
%         disp(['ip21: ',num2str(ip21)])
        disp(['----- | #P2M: ',num2str(ip21),' (out of: ',num2str(np2M),'), P2M: ',num2str(p21)]);
       
        for iMC=1:nMC
            
            state = zeros(NPop,1);
%             state = unifrnd(0,1,NPop,1)<rand;
            stArr = nan(NPop,tf);
            
            for t=1:tf
                stArr(:,t) = state;
                
%                 randArr = rand(NPop,1);
%                 posib20 = p20*state;
%                 posib21 = p21*not(state);
%                 
%                 state(randArr<posib20) = 0;
%                 state(randArr<posib21) = 1;
                
                randArr = rand(NPop,1);
                
                posib2msngr = p20*state;
                posib2explt = p21*not(state);
                
                expChange = randArr<posib2msngr;
                msgChange = randArr<posib2explt;
                
                state(expChange) = not(state(expChange));
                state(msgChange) = not(state(msgChange));
                
            end
            
            expltDur(iMC,:,ip20,ip21) = mean(stArr,2); % ratio of time each agent stayed in Exploiter state
            tmp = mean(stArr);
            numExpl(iMC,:,ip20,ip21) = tmp(1:nSkipSave:end);   % mean ratio of pop. doing exploitation
            tmp = diff(stArr,1,2);
            numEdges(iMC,:,ip20,ip21) = sum(abs(tmp),2);
            
        end
    end
end

fig = figure(1);
plot(stArr(1,:), 'LineWidth', 2.5)

fig.Position = [25.6667 198 1.8127e+03 420];
fig.Color = 'w';

ax = gca;
ax.YTick = [0 1];
ax.YTickLabel = {'Exploiter','Messenger'};
xlabel("Time Step")

set(gca,'FontSize',20)

nEdges = mean(numEdges(iMC,:,ip20,ip21),2);

figure(2)

plot([1:size(stArr,1)], stArr, 'LineWidth', 0.5)
disp("num edges: "+ nEdges)


