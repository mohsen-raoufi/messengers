clc; clear all;

% addpath('../../01_AES/scioi_01_aes/Functions/');
% addpath('../Functions/');
addpath('../');
addpath('../post_process/');

dateStr = date;
dateStr = replace(dateStr,'-','_');

nMC = 8;
tf = 50e3;
nTVars = 100;
nSkipSave = round(tf/nTVars);
NPop = 100;
linkThresh = 0.15; % the sensng range below which two nodes get connected
% linkThresh = 2*0.10; % the sensng range below which two nodes get connected
ArenaScale = 1.0;

% mainStrng = strcat('not_to_push_data/MCs_',dateStr,'__NPop_',num2str(NPop), '_Arena_',num2str(ArenaScale),'__tf_',num2str(tf/1000),'k__',...
%     'lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_',num2str(linkThresh),'__nMC_',num2str(nMC)); % _initENumMsngr

mainStrng = strcat('not_to_push_data/MCs_',dateStr,'__NPop_',num2str(NPop), '_Arena_',num2str(ArenaScale),'__tf_',num2str(tf/1000),'k__',...
    'lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_',num2str(linkThresh),'__nMC_',num2str(nMC),'_noNoise'); % _initENumMsngr
%% general config.s
neibZS = 5; % 0 for neigh: ZP, 1 for neigh: ZS, 2: ZPi+ZPj, 3: ZPi+ZSj, 4: ZSi+ZSj, 5: ZSi+ZPj


% landFunc = @(x,y) (x.^2 + y.^2);          % convex/bowl shape

% landFunc = @(x,y) (abs(x) + abs(y));      % pyramid shape
% land_func_name = 'pyramid';

landFunc = @(x,y) sqrt(x.^2 + y.^2);        % cone shape
land_func_name = 'cone';

% landFunc = @(x,y) (x.^2+y.^2).^0.25;      % square root of cone shape

% landFunc = @(x,y) x + y;
% land_func_name = 'ramp';

% landFunc = @(x,y) x;
% land_func_name = 'ramp_x';

% O = -[0 0.5];
% landFunc = @(x,y) ((1-(x-O(1))).^2)+(100*(((y-O(2))-((x-O(1)).^2)).^2));        % Rosenbrock Function from CEC 2005
% land_func_name = 'Rosenbrock';

% landFunc = @(x,y) sin(4*pi*x) + sin(4*pi*y);        % sine shape
% land_func_name = 'sine';
% landFunc = @(x,y) sin(3*pi*x) + sin(3*pi*y);        % sine_2 shape
% land_func_name = 'sine2';

% landFunc = @(x,y) sin(2*pi*x) + sin(2*pi*y);        % sine_2 shape
% land_func_name = 'sine2_2pix';

% % peakz = @(xx,yy) 3*(1-xx).^2.*exp(-(xx.^2) - (yy+1).^2) - 10*(xx/5 - xx.^3 - yy.^5).*exp(-xx.^2-yy.^2) - 1/3*exp(-(xx+1).^2 - yy.^2) ;
% % landFunc = @(x,y) peakz(2*x,2*y) ;        % peak function
% landFunc = @(x,y) peakz(x,y) ;        % peak function
% land_func_name = 'peakz_x_y';

% sig = 0.1; mu = 0.2;
% gausFun = @(r) 1/sig/sqrt(2*pi)*exp(-0.5*((r-mu)/sig).^2);
% landFunc = @(x,y) gausFun(sqrt(x.^2+y.^2));
% land_func_name = 'donut';

% k = 2*pi;
% landFunc = @(x,y) 2*((x) + (y)).^1 + sin(k*x) + sin(k*y);%

landBounds = ArenaScale*[-1 1; -1 1]; %
initBounds = landBounds;

x = linspace(landBounds(1,1),landBounds(1,2),100);
y = linspace(landBounds(2,1),landBounds(2,2),100);
[xx,yy] = meshgrid(x,y);
zz = landFunc(xx,yy);
% landFunc = @(x,y)landFunc(x,y);
landFunc = @(x,y) (landFunc(x,y)-min(zz(:)))/(max(zz(:))-min(zz(:)));
zz = landFunc(xx,yy);
zEnv = mean(zz(:));

%% parameters
Alpha = 0.99;
beta = 0.5;%0.01;   % forgetting factor on gradient: COLLECTIVE.
sigma = 0.001;
randStep0 = 1.0;
randStep2 = 0.001;
stepSize = 0.002;


% %% NEW PARAMETERS
alpha = 0.99;
beta = 0.9; % 0.5;   % forgetting factor on gradient: COLLECTIVE.
sigma = 0.0001; % 0.001
randStep0 = 1.0;
randStep2 = 0.25*0.0001;

% JUST FOR DEBUG
% p2expltArr = exp(-[1:20:50]./5);
% p2msngrArr = p2expltArr(1:end-1);


% p2expltArr = exp(-[1:2:50]./5);
% p2msngrArr = p2expltArr(1:end-1);

% p2expltArr = exp(-[1:4:100]./5);
% p2msngrArr = p2expltArr(1:end-1);

% Larger Parameter Set Regime
p2expltArr = exp(-[1:2:100]./5);
p2msngrArr = p2expltArr(1:end-1);

nVar1 = length(p2expltArr);
nVar2 = length(p2msngrArr);

strng = strcat(mainStrng,'_Params.mat');
save(strng);

%% init void variables
z1MeanArr  = nan(nVar1,nVar2,nTVars,nMC);
z1StdArr   = z1MeanArr;

posArr = nan(nVar1,nVar2,2,NPop,nTVars,nMC);
stateArr = nan(nVar1,nVar2,NPop,nTVars,nMC);
zpArr = nan(nVar1,nVar2,NPop,nTVars,nMC);
zsArr = zpArr;
% debugArr = stateArr;

wbar = waitbar(0, 'Starting');
for ip2E=1:nVar1 % for study on switching time
    % for iT2W=1:nVar1 % for study on switching time
    waitbar(ip2E/nVar1, wbar, sprintf('Progress: %d %%', floor(ip2E/nVar1*100)));
    p2explt = p2expltArr(ip2E);
    disp(['#P2E: ',num2str(ip2E),' (out of: ',num2str(nVar1),'), P2E: ',num2str(p2explt)]);

    %     time2waitInMsngr = time2waitInMsngrArr(iT2W);
    %     disp(['#t2Wait: ',num2str(iT2W),' (out of: ',num2str(nVar1),'), t2Wait: ',num2str(time2waitInMsngr)]);

    % iMC = 1;
    % parfor ip2M=1:nVar2
    
    for ip2M=1:nVar2
        p2msngr = p2msngrArr(ip2M);
        if rem(ip2M,10)==1
            disp(['----- | #P2M: ',num2str(ip2M),' (out of: ',num2str(nVar2),'), P2M: ',num2str(p2msngr)]);
        end

        parfor iMC=1:nMC
            

            pos = unifrnd(initBounds(1,1),initBounds(1,2),[1,NPop]);
            pos = [pos; unifrnd(initBounds(2,1),initBounds(2,2),[1,NPop])];

            [z1MeanArr(ip2E,ip2M,:,iMC), z1StdArr(ip2E,ip2M,:,iMC), posArr(ip2E,ip2M,:,:,:,iMC),...
                stateArr(ip2E,ip2M,:,:,iMC), zpArr(ip2E,ip2M,:,:,iMC), zsArr(ip2E,ip2M,:,:,iMC), ~] = ...
                funcEEM_Markov_new(landFunc,landBounds,pos,linkThresh,tf,nTVars,...
                neibZS,Alpha,beta,sigma,p2explt,p2msngr,randStep0,randStep2,stepSize);
            %
            %             [z1MeanArr(iT2W,ip2M,:,iMC), z1StdArr(iT2W,ip2M,:,iMC), posArr(iT2W,ip2M,:,:,:,iMC),...
            %                 stateArr(iT2W,ip2M,:,:,iMC), zpArr(iT2W,ip2M,:,:,iMC), zsArr(iT2W,ip2M,:,:,iMC), ~] = ...
            %                 funcEEM_constantWaitingTime(landFunc,landBounds,pos,linkThresh,tf,nTVars,...
            %                 neibZS,Alpha,beta,sigma,time2waitInMsngr,p2msngr,randStep0,randStep2,stepSize);
        end
    end
end
close(wbar);
close all;


strng = mainStrng;
save(strng);

steadyStateFileMaker();

beep