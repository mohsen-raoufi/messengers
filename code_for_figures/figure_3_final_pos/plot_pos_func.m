function plot_pos_func(pos, z, state, xx, yy, zz, zEnv, NPop, ArenaScale, linkThresh)

NRobots = NPop;

plot(nan,nan)

set(gca,'fontsize',12);
set(gcf,'color','w');

hold on;

% Red to purple to blue
% nColors = 200;
% custom_cmap_map = [linspace(0,1,nColors); zeros(1,nColors); linspace(1,0,nColors)]'.^0.5;
% cMAP = colormap(custom_cmap_map);

% abyss and copper merged
% nColors = 200;
% cmap_abyss = flipud(abyss(nColors/2));
% cmap_copper = copper(nColors/2);
% % custom_cmap_map = [cmap_abyss.^0.5; cmap_copper.^0.5];
% cmap_abyss = cmap_abyss./cmap_abyss(1,:);
% cmap_copper = cmap_copper./cmap_copper(end,:);
% custom_cmap_map = [cmap_abyss; cmap_copper];
% cMAP = colormap(custom_cmap_map);

% load from external file
custom_cmap_map = load('Berlin_RdBu_denser.mat');
custom_cmap_map = custom_cmap_map.RdBu;
cMAP = colormap(custom_cmap_map);

% cMAP = abyss(100);

maxZ = max(zz(:));
minZ = min(zz(:));
lenZ = maxZ - minZ;

ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
set(gca,'FontSize',14)

cbar.Limits = [0, 1];
caxis([0, 1]);

rng('shuffle');


% plot the precision error as a shade!
radii = sqrt(pos(1,:).^2 + pos(2,:).^2);
mean_rad = mean(radii);
bnd_w = std(radii);
hi_bnd = mean_rad + bnd_w;
low_bnd = max(mean_rad - bnd_w,0.001);
circle_0 = nsidedpoly(1000,'Center',[0,0],'Radius',mean_rad);
circle_1 = nsidedpoly(1000,'Center',[0,0],'Radius',hi_bnd);
circle_2 = nsidedpoly(1000,'Center',[0,0],'Radius',low_bnd);

p0 = plot(circle_0,'FaceColor','w','EdgeColor','b');
p1 = plot(circle_1,'FaceColor','k','EdgeColor','w');
p2 = plot(circle_2,'FaceColor','w','EdgeColor','w');
alpha(p1,0.1);
alpha(p2,1);


axis square
axis([-1 1 -1 1]*ArenaScale)
% title(['Time Step: ',num2str(time)]);
contour(xx,yy,zz,[zEnv zEnv],'-.k','LineWidth',1.5,'EdgeAlpha',0.5)

% time = size(posArr,5);
%
% tmp = posArr(1,1,:,:,time,1);
% pos = reshape(tmp,[2,size(tmp,4)]);
%
% tmp = zpArr(1,1,:,time,1);
% z = reshape(tmp,[1,size(tmp,3)]);
%
% tmp = stateArr(1,1,:,time,1);
% state = reshape(tmp,[1,size(tmp,3)]);

time = 1;
[sc, sc2] = plot_pos_net_0(pos, z, state, cMAP, minZ, lenZ, NRobots, linkThresh, true, 70);

axis off
box on
% box off
xlm = xlim;
ylm = ylim;
line(xlm,[ylm(2), ylm(2)],'Color','k')
line(xlm,[ylm(1), ylm(1)],'Color','k')
line([xlm(2), xlm(2)],ylm,'Color','k')
line([xlm(1), xlm(1)],ylm,'Color','k')

contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2)
colormap(cMAP)

end




