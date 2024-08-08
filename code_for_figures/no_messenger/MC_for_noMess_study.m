clc; 
close all;
clear all;

% addpath('../');
addpath('../functions/')
% parpool(20);

dateStr = date;
dateStr = replace(dateStr,'-','_');

nMC = 40;
tf = 50e3;
nTVars = 20;
nSkipSave = round(tf/nTVars);
NPop = 100;

linkThreshArr = [0:0.005:0.5, 0.51:0.01:1]; %linspace(0,2,100);

ArenaScale = 1.0;

mainStrng = strcat('data/MCs_',dateStr,'__NPop_',num2str(NPop), '_Arena_',num2str(ArenaScale),'__tf_',num2str(tf/1000),'k__',...
    'cone__BasicMarkov__NoMessenger_sensRang_Study__nMC_',num2str(nMC)); % _initENumMsngr
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
% randStep2 = 0.04;
stepSize = 0.002;
% stepSize = 0.001;

% %% NEW PARAMETERS
alpha = 0.99;
beta = 0.9; % 0.5;   % forgetting factor on gradient: COLLECTIVE.
sigma = 0.0001; % 0.001
randStep0 = 1.0;
randStep2 = 0.25*0.0001;


p2explt = 1.0;
p2msngr= 0;


var2Study = linkThreshArr;
nVar = length(var2Study);

strng = strcat(mainStrng,'_Params.mat');
save(strng);

%% init void variables
z1MeanArr  = nan(nVar,nTVars,nMC);
z1StdArr   = z1MeanArr;

posArr = nan(nVar,2,NPop,nTVars,nMC);
stateArr = nan(nVar,NPop,nTVars,nMC);
zpArr = nan(nVar,NPop,nTVars,nMC);
zsArr = zpArr;
posArr_init = nan(nVar,2,NPop,1,nMC);

wbar = waitbar(0, 'Starting');
for iVar=1:nVar % for study on var to study
    waitbar(iVar/nVar, wbar, sprintf('Progress: %d %%', floor(iVar/nVar*100)));
    
    linkThresh = linkThreshArr(iVar);
    parfor iMC=1:nMC

        pos = unifrnd(initBounds(1,1),initBounds(1,2),[1,NPop]);
        pos = [pos; unifrnd(initBounds(2,1),initBounds(2,2),[1,NPop])];
        
        posArr_init(iVar,:,:,1,iMC) = pos;

        [z1MeanArr(iVar,:,iMC), z1StdArr(iVar,:,iMC), posArr(iVar,:,:,:,iMC),...
            stateArr(iVar,:,:,iMC), zpArr(iVar,:,:,iMC), zsArr(iVar,:,:,iMC), ~] = ...
            funcEEM_Markov_new(landFunc,landBounds,pos,linkThresh,tf,nTVars,...
            neibZS,Alpha,beta,sigma,p2explt,p2msngr,randStep0,randStep2,stepSize);
    end
end
close(wbar);
close all;

strng = mainStrng;
save(strng);

% plot_postPro_for_linkThresh();

beep