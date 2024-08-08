% funcEEM_Markov(landFunc,LB,pos,linkThresh,tf,nTVars,neibZS,alpha,beta,sigma,p2explt,p2msngr,randStep0,randStep2,stepSize)
clc; clear all; 
close all

bool_plot_online = false;

tf = 2.0e3;
nTVars = tf/1;
nSkipSave = round(tf/nTVars);
% NPop = 100;
NPop = 200;


alpha = 0.99;
beta = 0.8; % 0.5;   % forgetting factor on gradient: COLLECTIVE.
sigma = 0.002; % 0.001
randStep0 = 1.0;
randStep2 = 0.02;
% stepSize = 0.002;
stepSize = 0.0015;

p2explt = 1.0; 0.01;
p2msngr = 0.0; 0.91;

% linkThresh = 0.15; % the sensng range below which two nodes get connected
linkThresh = 0.40; % the sensng range below which two nodes get connected
ArenaScale = 0.7; 1.0;
landBounds = ArenaScale*[-1 1; -1 1]; %
initBounds = landBounds;
x = linspace(landBounds(1,1),landBounds(1,2),100);
y = linspace(landBounds(2,1),landBounds(2,2),100);
[xx,yy] = meshgrid(x,y);
LB = landBounds;

% landFunc = @(x,y) (x.^2 + y.^2);          % convex/bowl shape

% landFunc = @(x,y) (abs(x) + abs(y));      % pyramid shape
% land_func_name = 'pyramid';

% landFunc = @(x,y) sqrt(x.^2 + y.^2);        % cone shape
% land_func_name = 'cone';

% landFunc = @(x,y) (x.^2+y.^2).^0.25;      % square root of cone shape
% land_func_name = 'cone_sqrt';

% landFunc = @(x,y) x + y;
% land_func_name = 'ramp';

landFunc = @(x,y) x;
land_func_name = 'ramp_x';

% O = -[0 0.5];
% landFunc = @(x,y) ((1-(x-O(1))).^2)+(100*(((y-O(2))-((x-O(1)).^2)).^2));        % Rosenbrock Function from CEC 2005
% land_func_name = 'Rosenbrock';

% landFunc = @(x,y) sin(4*pi*x) + sin(4*pi*y);        % sine shape
% land_func_name = 'sine';
% landFunc = @(x,y) sin(3*pi*x) + sin(3*pi*y);        % sine_2 shape
% land_func_name = 'sine2';

% landFunc = @(x,y) sin(2*pi*x) + sin(2*pi*y);        % sine_2 shape
% land_func_name = 'sine2_2pix';

% peakz = @(xx,yy) 3*(1-xx).^2.*exp(-(xx.^2) - (yy+1).^2) - 10*(xx/5 - xx.^3 - yy.^5).*exp(-xx.^2-yy.^2) - 1/3*exp(-(xx+1).^2 - yy.^2) ;
% % landFunc = @(x,y) peakz(2*x,2*y) ;        % peak function
% landFunc = @(x,y) peakz(x,y) ;        % peak function
% land_func_name = 'peakz_x_y';

% sig = 0.1; mu = 0.2;
% gausFun = @(r) 1/sig/sqrt(2*pi)*exp(-0.5*((r-mu)/sig).^2);
% landFunc = @(x,y) gausFun(sqrt(x.^2+y.^2));
% land_func_name = 'donut';

% k = 2*pi;
% landFunc = @(x,y) 2*((x) + (y)).^1 + sin(k*x) + sin(k*y);%

zz = landFunc(xx,yy);

% landFunc = @(x,y)landFunc(x,y);
landFunc = @(x,y) (landFunc(x,y)-min(zz(:)))/(max(zz(:))-min(zz(:)));
zz = landFunc(xx,yy);
zEnv = mean(zz(:));

pos = unifrnd(initBounds(1,1),initBounds(1,2),[1,NPop]);
pos = [pos; unifrnd(initBounds(2,1),initBounds(2,2),[1,NPop])];

% swch = ones(1,NPop); % 0: Messenger, 1: Exploiter
ExpNExploiter = p2explt/(p2explt+p2msngr);
tmp = rand(1,NPop);
swch = tmp<ExpNExploiter;

randStep = randStep0*ones(1,NPop);

neibZS = 5; % 0 for neigh: ZP, 1 for neigh: ZS, 2: ZPi+ZPj, 3: ZPi+ZSj, 4: ZSi+ZSj, 5: ZSi+ZPj


%% Initialization
[~,ind] = sort(pos(2,:));
pos = pos(:,ind);

Adjc = zeros(NPop); % Adjacancy matrix, previously called Links!
for jj=1:NPop-1
    for ii=jj+1:NPop
        if(norm(pos(:,jj)-pos(:,ii))<linkThresh)
            Adjc(ii,jj) = 1;
        end
    end
end
Adjc = Adjc + Adjc';

cumulated_Adjc = Adjc;

ZR = landFunc(pos(1,:),pos(2,:));
ZS = ZR + sigma*randn(1,NPop);
ZP = ZS;
%% Parameters
% zz = rescale(zz);
Orient = unifrnd(-pi,pi,[1,NPop]);
rndOrient = unifrnd(-pi,pi,[1,NPop]);

% Walk = zeros(2,NPop);
% grad = zeros(2,NPop);
Walk = 0.0001*unifrnd(-1,1,2,NPop);
grad = 0.0001*unifrnd(-1,1,2,NPop);
prevgrad = grad;
curGrad = grad;

objFArrOld = 0;

% %% Variable to return
nSkipSave = round(tf/nTVars);
zpArr = nan(NPop,nTVars);
zsArr = zpArr;
z1Arr = nan(1,nTVars);
z1StdArr = z1Arr;
stArr = nan(NPop,nTVars);
posArr = nan(2,NPop,nTVars);
debugArr = stArr;

degr = sum(Adjc);
saveCtr = 0;
rng('shuffle');

if(bool_plot_online)
    fig = figure(1);
    hold on;
    sct = scatter(pos(1,:), pos(2,:), 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0);
    axis([-1 1 -1 1])
    axis square
end
%% Main loop
for time=1:tf
    %     swch = 2*(time>tSwitch)-1;

    %     stArr(:,saveCtr) = swch;

    prevSwch = swch;



    ZR = landFunc(pos(1,:),pos(2,:));
    %     zz = landFunc(xx-Oo(1),yy-Oo(2));

    ZS = ZR + sigma*randn(1,NPop); % << Adding Noise

    % % % % %     if(time==tSwitch)
    % % % % %         swch = ones(size(swch));
    % % % % %         ZP = ZS;
    % % % % %     end

    %     if(rem(time,50)==0)
    randArr = rand(1,NPop);

    posib2msngr = p2msngr*swch;         % Basic one, constant
    posib2explt = p2explt*not(swch);

    %         ts = time/tf;
    %         posib2msngr = p2msngr*swch*(1-ts); % temporally decaying prob.
    %         posib2explt = p2explt*not(swch);

    expChange = randArr<posib2msngr;
    msgChange = randArr<posib2explt;


    %         swch(expChange) = and((degr(swch(expChange))>1), not(swch(expChange)) );
    %         swch(msgChange) = not(swch(msgChange));



    swch(expChange) = not(swch(expChange));
    swch(msgChange) = not(swch(msgChange));


    %         swch(randArr<posib2msngr) = 0;
    %         swch(randArr<posib2explt) = 1;

    %         swch = zeros(1,NPop);
    
    %     end

    %     if(time>1000)
    %         aaa = 2;
    %
    %     end

    randStep(swch==1) = randStep2;
    randStep(swch==0) = randStep0;

    %     timeConst = 100;
    %     if(rem(time,timeConst)==1)
    %         rFlip = 0.0;
    %         nFlip = ceil(rFlip*NRobots);
    %         indSwFlip = randperm(NRobots,nFlip);
    %         swch(indSwFlip) = not(swch(indSwFlip));
    %         randStep(indSwFlip) = randStep0;
    %     elseif(rem(time,timeConst)==0)
    %         swch(indSwFlip) = not(swch(indSwFlip));
    %         randStep(indSwFlip) = randStep2;
    %     end


    %     swch = time>tSwitch;

    justSwitched = swch~=prevSwch;
    grad(:,justSwitched) = rand(2,sum(justSwitched));
    prevgrad(:,justSwitched) = rand(2,sum(justSwitched));

    z1 = mean(ZS,'all');
    stdZ1 = std(ZS);



    %% Gradient Search
    degr = sum(Adjc);
    switch neibZS
        case 0
            neigh = ZP*(Adjc)./degr;
        case 1
            neigh = ZS*(Adjc)./degr; % ZP
        case 2
            neigh = (ZP+ZP*(Adjc))./(degr+1);
        case 3
            neigh = (ZP+ZS*(Adjc))./(degr+1);
        case 4
            neigh = (ZS+ZS*(Adjc))./(degr+1);
        case 5
            neigh = (ZS+ZP*(Adjc))./(degr+1);
    end

    %     switch indMem
    %         case 0
    % %             ZP = 0.01*ZP + 0.99*ZS;
    %             ZP = 0.99*ZP + 0.01*ZS;
    %         case 1
    %             minZ = min(minZ,ZS);
    %             maxZ = max(maxZ,ZS);
    %             ZminMax = 0.5*(minZ+maxZ);
    % %             ZP = 0.5*(minZ+maxZ);
    % %             min max
    %         case 2
    %             ZP = (time*ZP+ZS)./(time+1);
    %         case 3
    %         case 4
    %     end


    %     ZP(swch~=1) = beta_1*ZminMax(swch~=1) + (1-beta_1)*neigh(swch~=1);
    %     ZP(swch~=1) = ZminMax(swch~=1);
    %     ZP(swch==1) = beta_1*ZP(swch==1) + (1-beta_1)*neigh(swch==1);
    %     ZP_ref = beta_2*ZP_ref + (1-beta_2)*ZP;

    %     ZP(swch~=1) = beta_1*ZminMax(swch~=1) + (1-beta_1)*neigh(swch~=1);
    %     ZP(swch~=1) = ZminMax(swch~=1);

    ZP(swch~=1) = ZP(swch~=1);                                          % non-exploiting agents = messengers
    ZP(swch==1) = alpha*ZP(swch==1) + (1-alpha)*neigh(swch==1);       % exploiting agents

    objFArr = ZS;                                                       % For extended pseudo-gradient
    %     objFArr = swch.*abs(ZS-ZP);                                         % for the basic pseudo-gradient

    diffObjF = objFArr - objFArrOld;
    objFArrOld = objFArr;

    tmpAng = atan2(Walk(2,:),Walk(1,:));
    

    % curGrad = [diffObjF; diffObjF]./Walk;                         % Avoid using this. It will lead to some weird diagonal-only motions due to numerical issues
        
    curGrad = diffObjF.*[cos(tmpAng); sin(tmpAng)];

    grad = beta*prevgrad + (1-beta)*curGrad;

    prevgrad = curGrad;

    % % %     if(time>tSwitch)
    % % %         grad = Gz;
    % % %     else
    % %     if(time<tSwitch)
    % %         grad=rand(2,NPop);
    % %         prevgrad=rand(2,NPop);
    % %     end

    if(sum(isnan(grad),'all'))
        disp('grad is NAN');
        grad(isnan(grad)) = 0;
        if(time>2), disp('**'); end
    end

    gradNormed = normc(grad);

    % %% Two different random motions for the two states
    rndOrient(swch~=1) = rndOrient(swch~=1) + 0.75*0.1*unifrnd(-pi,pi,[1,sum(swch~=1)]);
    %         rndOrient(swch~=1) = unifrnd(-pi,pi,[1,sum(swch~=1)]);
    rndOrient(swch==1) = unifrnd(-pi,pi,[1,sum(swch==1)]);

    rndOrient = rndOrient + 0.1*unifrnd(-pi,pi,[1,NPop]);


    rndWalkU = [cos(rndOrient); sin(rndOrient)];


    if(time>2)

                % Walk = ((1-randStep).*grad);
                % Walk = Walk + (randStep.*rndWalkU).*[vecnorm(Walk); vecnorm(Walk)];
                % Walk = -stepSize*normc(Walk);


        % %         % Original 1st Gradient Descent
                % Walk = ((1-randStep).*gradNormed);
                % Walk = Walk + (randStep.*rndWalkU);
                % Walk = -stepSize*normc(Walk);

        Walk = ((1-randStep).*(ZS-neigh).*(degr./(1+degr)).*(grad));
        Walk = Walk + randStep.*rndWalkU.*([vecnorm(Walk); vecnorm(Walk)]);
        Walk = -stepSize.*normc(Walk);

    else
        Walk = stepSize*rndWalkU;
    end

    Walk(:,swch==0) = stepSize.*rndWalkU(:,swch==0);


    if(sum(isnan(Walk),'all'))
        indxNaN = find(isnan(Walk));
        Walk(indxNaN) = stepSize*randStep*rndWalkU(indxNaN);
        fprintf("Walk is NAN! Number of robots: %i, time: %i\n",length(indxNaN),time);
    end

    if(sum(sum(Adjc)==0)) % lonely agents just walk randomly
        indxDetach = find(sum(Adjc)==0);
        Walk(:,indxDetach) = stepSize*rndWalkU(:,indxDetach);
    end

    %     Walk = min(max(Walk,-0.01),0.01); % <<<<<< Never use min max on Walk
    nWalk = vecnorm(Walk);
    nWalkN = min(nWalk,stepSize);
    Walk = Walk./nWalk.*nWalkN;

    prevPos = pos;
    pos = pos + Walk;

    if(isnan(grad) + isnan(Walk))
        disp('HOOY');
    end

    pos(1,:) = min(max(pos(1,:),LB(1,1)),LB(1,2));
    pos(2,:) = min(max(pos(2,:),LB(2,1)),LB(2,2));


    % %% Wall reflecting boundary condition
    punishMargin = 0.05;
    if( (min(pos(1,:))<(LB(1,1)+punishMargin*linkThresh)) || (max(pos(1,:))>(LB(1,2)-punishMargin*linkThresh)) ) % x-axis
        indPassed = find(pos(1,:)<(LB(1,1)+punishMargin*linkThresh));
        pos(1,indPassed) = pos(1,indPassed) + 2*abs(Walk(1,indPassed));%- 2*Walk(1,nPassed);
        rndrnd = unifrnd(-pi,pi,[NPop,1]);
        rndOrient(indPassed) = rndrnd(indPassed);%rndOrient(nPassed) + 0.3;

        indPassed = find(pos(1,:)>(LB(1,2)-punishMargin*linkThresh));
        pos(1,indPassed) = pos(1,indPassed) - 2*abs(Walk(1,indPassed));%- 2*Walk(1,nPassed);
        rndOrient(indPassed) = rndrnd(indPassed);%rndOrient(nPassed) + 0.3;

        Walk(:,indPassed) = stepSize*[cos(rndOrient(indPassed)); sin(rndOrient(indPassed))];
        %         disp('punished');
    end

    if( (min(pos(2,:))<(LB(2,1)+punishMargin*linkThresh)) || (max(pos(2,:))>(LB(2,2)-punishMargin*linkThresh)) ) % y-axis
        indPassed = find(pos(2,:)<(LB(2,1)+punishMargin*linkThresh));
        pos(2,indPassed) = pos(2,indPassed) + 2*abs(Walk(2,indPassed));%- 2*Walk(2,nPassed);
        %         rndOrient(nPassed) = rndOrient(nPassed) + 0.3;
        rndrnd = unifrnd(-pi,pi,[NPop,1]);
        rndOrient(indPassed) = rndrnd(indPassed);%rndOrient(nPassed) + 0.3;

        indPassed = find(pos(2,:)>(LB(2,2)-punishMargin*linkThresh));
        pos(2,indPassed) = pos(2,indPassed) - 2*abs(Walk(2,indPassed));%- 2*Walk(2,nPassed);
        %         rndOrient(nPassed) = rndOrient(nPassed) + 0.3;
        rndOrient(indPassed) = rndrnd(indPassed);
        Walk(:,indPassed) = stepSize*[cos(rndOrient(indPassed)); sin(rndOrient(indPassed))];
        %         disp('punished 2');
    end

    Adjc = zeros(NPop);
    for mm=1:NPop-1
        for nn=mm+1:NPop
            if(norm(pos(:,mm)-pos(:,nn))<linkThresh)
                Adjc(nn,mm) = 1;
            end
        end
    end
    Adjc = Adjc + Adjc';

    cumulated_Adjc = and(Adjc,cumulated_Adjc);

    if(bool_plot_online)

        if(rem(time,10)==0)
            % delete(sct)
            % clr = (time/tf).^0.5*[1 0 0]; % red
            clr = (time/tf).^0.5*[0 0.5 1]; % blue
            sct = scatter(pos(1,:), pos(2,:), 5, clr, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0);
            %     scatter(posArr(1,:,t), posArr(2,:,t), 10, clr, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0);
            drawnow();
        end

    end

    if rem(time,nSkipSave)==0
        saveCtr = saveCtr + 1;

        posArr(:,:,saveCtr) = pos;

        zpArr(:,saveCtr) = ZP;
        zsArr(:,saveCtr) = ZS;
        z1Arr(saveCtr) = z1;
        z1StdArr(saveCtr) = stdZ1;

        stArr(:,saveCtr) = swch;

        debugArr(:,saveCtr) = tmpAng;

    end

end

disp("finished simulation!!");

if (~bool_plot_online)

    fig = figure(2);
    fig.Position = [680 458 600 600];
    hold on;
    axis(ArenaScale*[-1 1 -1 1])
    % axis([-0 2 -0 2])
    axis square



    for time=1:nTVars
        tmp = posArr(:,:,time);
        pos = reshape(tmp,[2,NPop]);

        % delete(sct)
        % clr = (time/nTVars).^1*[1 0 0];
        clr = (time/nTVars).^0.5*[0 0.5 1]; % blue
        sct = scatter(pos(1,:), pos(2,:), 7, clr, 'filled', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0);
        %     scatter(posArr(1,:,t), posArr(2,:,t), 10, clr, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0);
        % drawnow();


    end


    set(gca,'Color','w')
    set(gcf,'Color','w')
    axis off
end

save_str = land_func_name+"_noisy_N_"+ num2str(NPop);
savefig(save_str + "_fig.fig")
save(save_str + "_data.mat")
exportgraphics(fig,save_str+"_trace.png",'resolution',300)

