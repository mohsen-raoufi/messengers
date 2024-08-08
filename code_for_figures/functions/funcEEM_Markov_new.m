function [z1Arr, z1StdArr, posArr, stArr, zpArr, zsArr, debugArr] = ...
    funcEEM_Markov_new(landFunc,LB,pos,linkThresh,tf,nTVars,neibZS,alpha,beta,sigma,p2explt,p2msngr,randStep0,randStep2,stepSize)

NPop = size(pos,2);
% swch = ones(1,NPop); % 0: Messenger, 1: Exploiter
ExpNExploiter = p2explt/(p2explt+p2msngr);
tmp = rand(1,NPop);
swch = tmp<ExpNExploiter;

randStep = randStep0*ones(1,NPop);

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
grad = unifrnd(-1,1,2,NPop);
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
%% Main loop
for time=1:tf
    %     swch = 2*(time>tSwitch)-1;
    if rem(time,nSkipSave)==1
        saveCtr = saveCtr + 1;
    end
    %     stArr(:,saveCtr) = swch;

    prevSwch = swch;
    posArr(:,:,saveCtr) = pos;


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
    stArr(:,saveCtr) = swch;
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
    grad(:,justSwitched) = unifrnd(-1,1,2,sum(justSwitched)); % newer version
    % grad(:,justSwitched) = rand(2,sum(justSwitched)); % old version
    prevgrad(:,justSwitched) = unifrnd(-1,1,2,sum(justSwitched)); % newer version
    % prevgrad(:,justSwitched) = rand(2,sum(justSwitched)); % old version

    z1 = mean(ZS,'all');
    stdZ1 = std(ZS);

    zpArr(:,saveCtr) = ZP;
    zsArr(:,saveCtr) = ZS;
    z1Arr(saveCtr) = z1;
    z1StdArr(saveCtr) = stdZ1;

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
    debugArr(:,saveCtr) = tmpAng;

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
    % rndOrient(swch~=1) = unifrnd(-pi,pi,[1,sum(swch~=1)]);

    rndOrient(swch==1) = unifrnd(-pi,pi,[1,sum(swch==1)]);

    %         rndOrient = rndOrient + 0.1*0.1*unifrnd(-pi,pi,[1,NPop]);


    rndWalkU = [cos(rndOrient); sin(rndOrient)];


    if(time>2)

        %         Walk = ((1-randStep).*grad);
        %         Walk = Walk + (randStep.*rndWalkU).*[vecnorm(Walk); vecnorm(Walk)];
        %         Walk = -stepSize*normc(Walk);


        % %         % Original 1st Gradient Descent
        %         Walk = ((1-randStep).*gradNormed);
        %         Walk = Walk + (randStep.*rndWalkU);
        %         Walk = -stepSize*normc(Walk);

        Walk = ((1-randStep).*(ZS-neigh).*(degr./(1+degr)).*(grad));
        Walk = Walk + (randStep.*rndWalkU);
        % Walk = Walk + randStep.*rndWalkU.*([vecnorm(Walk); vecnorm(Walk)]);
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
        %
        % disp('punished 2');
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


    %     if(rem(time,1000)==0)
    %         disp(time);
    %     end

end

zpArr(:,end) = ZP;
zsArr(:,end) = ZS;
z1Arr(end) = z1;
z1StdArr(end) = stdZ1;
posArr(:,:,end) = pos;
stArr(:,end) = swch;
