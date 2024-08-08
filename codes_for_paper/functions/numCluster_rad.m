function groupIndxArr = numCluster_rad(xPos,yPos,rad)

% Poses = Arr;
NROBOTS = length(xPos);
Arr = [xPos,yPos,zeros(NROBOTS,1)];
% scatter(Arr(:,1),Arr(:,2));

for numel_i=1:size(Arr,1)
    pos_i = Arr(numel_i,1:2);
    if(Arr(numel_i,3)==0)
        groupIndex = max(Arr(:,3))+1;
        Arr(numel_i,3) = groupIndex;
    else
        groupIndex = Arr(numel_i,3);
    end
%     Arr(1:13,3)'
%     aa=2;
    for numel_j=numel_i+1:size(Arr,1)
%         if(numel_i==6)
%             if(numel_j==13)
%                 aa = 0; end
%             aa = 2;
%         end
        pos_j = Arr(numel_j,1:2);
        r_ij =  pos_i - pos_j;
        isSensible = isInSenseRange(r_ij,rad);
        if(isSensible)
            if (Arr(numel_j,3)==0)
                Arr(numel_j,3) = groupIndex;
            else
                repIndx = Arr(:,3)==groupIndex;
                groupIndex = Arr(numel_j,3);
                Arr(repIndx,3) = groupIndex;
            end
        end
    end
    
end

groupIndxArr = Arr(:,3);
uniqG = unique(groupIndxArr);
for jj=1:length(uniqG)
    repIndx = groupIndxArr (:) == uniqG(jj);
    groupIndxArr(repIndx) = jj;
end
% figure(1)
% plot(nan,nan);
% hold on;
% Arr(:,3)
% cmap = jet(max(Arr(:,3)));
% [~, ~, group_idx] = unique(Arr(:,3));
% col_for_val = cmap(group_idx(:), :);
% scatter(Arr(:,1),Arr(:,2),[],col_for_val,'filled')


% for jj=1:max(Arr(:,3))
%     indx = find(Arr(:,3)==jj);
%     k = boundary(Arr(indx,1),Arr(indx,2));
%     plot(Arr(indx(k),1),Arr(indx(k),2),'LineStyle','-','Color',col_for_val(indx(1),:));
% end

end