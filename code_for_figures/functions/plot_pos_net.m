function [sc, sc2] = plot_pos_net(posArr, zpArr, stateArr, time, cMAP, minZ, lenZ, NRobots, linkThresh, show_network, marker_size, bool_from_arr)
if(nargin<12 || ~bool_from_arr)

    tmp = posArr(1,1,:,:,time,1);
    pos = reshape(tmp,[2,size(tmp,4)]);

    tmp = zpArr(1,1,:,time,1);
    z = reshape(tmp,[1,size(tmp,3)]);

    tmp = stateArr(1,1,:,time,1);
    state = reshape(tmp,[1,size(tmp,3)]);
else
    pos = posArr;
    z = zpArr;
    state = stateArr;
end
indCol = round((z-minZ)/(lenZ)*(length(cMAP)-1));
indCol = max(indCol,1);
% disp(indCol)
col = cMAP(indCol,:);



if(show_network)

    Adjc = zeros(NRobots); % Adjacancy matrix, previously called Links!
    for jj=1:NRobots-1
        for ii=jj+1:NRobots
            if(norm(pos(:,jj)-pos(:,ii))<linkThresh)
                Adjc(ii,jj) = 1;
            end
        end
    end
    Adjc = Adjc + Adjc';

    for jj=1:NRobots
        for ii=jj+1:NRobots
            if Adjc(jj,ii) == 1
                plot([pos(1,jj) pos(1,ii)],[pos(2,jj) pos(2,ii)],'Color',[153 255 204 255]/255,'LineWidth',0.75+0.25);
            end
        end
    end

end

sc = scatter(pos(1,state==1),pos(2,state==1),marker_size,col(state==1,:),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
%     sc2 = scatter(pos(1,state==0),pos(2,state==0),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
if(show_network)
    sc2 = scatter(pos(1,state==0),pos(2,state==0),1.5*marker_size,col(state==0,:),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor','r', 'LineWidth',1.5);
else
    sc2 = scatter(pos(1,state==0),pos(2,state==0),marker_size,col(state==0,:),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
end
end