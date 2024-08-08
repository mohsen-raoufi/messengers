function strng = run_one_config_no_msngr_func(linkThresh)
clc;


i_p2e = 1;
i_p2m = 49;

% No Messenger [1, 49], [21, 30] , [36,45]

% addpath('../../01_AES/scioi_01_aes/Functions/');
addpath('../functions/');
% addpath('../');
% addpath('../../post_process/');

dateStr = date;
dateStr = replace(dateStr,'-','_');

nMC = 1;
tf = 50e3;
nTVars = 10; tf;
% tf = 50e3;
% nTVars = 1000; tf;
nSkipSave = round(tf/nTVars);
NPop = 100;
% linkThresh = 0.15; % the sensng range below which two nodes get connected
% linkThresh = 2*0.10; % the sensng range below which two nodes get connected
ArenaScale = 1.0;

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


% Alpha = 0.99;
% beta = 0.5;%0.01;   % forgetting factor on gradient: COLLECTIVE.
% sigma = 0.001;
% randStep0 = 1.0;
% randStep2 = 0.001;
% stepSize = 0.002;


% Alpha = 0.99;
% beta = 0.8;%0.01;   % forgetting factor on gradient: COLLECTIVE.
% sigma = 0.01;
% randStep0 = 1.0;
% randStep2 = 0.2*0.05;
% stepSize = 0.001;

% JUST FOR DEBUG
% p2expltArr = exp(-[1:20:50]./5);
% p2msngrArr = p2expltArr(1:end-1);


% p2expltArr = exp(-[1:2:50]./5);
% p2msngrArr = p2expltArr(1:end-1);

% Larger Parameter Set Regime
p2expltArr_large = exp(-[1:2:100]./5);
p2msngrArr_large = p2expltArr_large(1:end-1);



% 6,7 --- 9, 10

% Single parameter check
p2expltArr = p2expltArr_large;
p2msngrArr = p2msngrArr_large;
p2explt = p2expltArr(i_p2e);
p2msngr = p2msngrArr(i_p2m);

nVar1 = 1;
nVar2 = 1;


mainStrng = strcat('data/single_run_',dateStr,'__NPop_',num2str(NPop), '_Arena_',num2str(ArenaScale),'__tf_',num2str(tf/1000),'k__',...
    'lowRand4Explt__',land_func_name,'__BasicMarkov__initENumMsngr_sensRang_',num2str(linkThresh) ...
    ,'__i_p2e_',num2str(i_p2e), '__i_p2m_',num2str(i_p2m),'.mat'); % _initENumMsngr

strng = strcat(mainStrng(1:end-4),'_Params.mat');
save(strng);

%% init void variables
z1MeanArr  = nan(nVar1,nVar2,nTVars,nMC);
z1StdArr   = z1MeanArr;

posArr = nan(nVar1,nVar2,2,NPop,nTVars,nMC);
stateArr = nan(nVar1,nVar2,NPop,nTVars,nMC);
zpArr = nan(nVar1,nVar2,NPop,nTVars,nMC);
zsArr = zpArr;
% debugArr = stateArr;




pos = unifrnd(initBounds(1,1),initBounds(1,2),[1,NPop]);
pos = [pos; unifrnd(initBounds(2,1),initBounds(2,2),[1,NPop])];

iMC = 1;

[z1MeanArr(1,1,:,iMC), z1StdArr(1,1,:,iMC), posArr(1,1,:,:,:,iMC),...
    stateArr(1,1,:,:,iMC), zpArr(1,1,:,:,iMC), zsArr(1,1,:,:,iMC), ~] = ...
    funcEEM_Markov_new(landFunc,landBounds,pos,linkThresh,tf,nTVars,...
    neibZS,Alpha,beta,sigma,p2explt,p2msngr,randStep0,randStep2,stepSize);

close all;


strng = mainStrng;
save(strng, '-v7.3');

% figure(101);
% plot_trace_low_graphics()
% 
% figure(1);
% res_video_maker_paper();


% steadyStateFileMaker();
% 
% figure(888);
% plot(nan,nan);
% hold on
% for t=1:nTVars
%     plot(reshape(stateArr(1,1,:,t),[],1),'Color',[t/nTVars,0,1-t/nTVars])
% end
%     hold on
%     plot(reshape(stateArr(1,1,:,end),[],1),'b')

beep