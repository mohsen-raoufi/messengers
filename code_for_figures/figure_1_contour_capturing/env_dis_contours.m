clc;
% clear all;
close all;

fig = figure(3);
fig.Position = [680 458 500 500];
% landFunc = @(x,y) (x.^2 + y.^2);          % convex/bowl shape

landFunc = @(x,y) (abs(x) + abs(y));      % pyramid shape
land_func_name = 'pyramid';

landFunc = @(x,y) sqrt(x.^2 + y.^2);        % cone shape
land_func_name = 'cone';

% landFunc = @(x,y) (x.^2+y.^2).^0.25;      % square root of cone shape

% landFunc = @(x,y) x + y;
% land_func_name = 'ramp';
% 
landFunc = @(x,y) x;
land_func_name = 'ramp_x';
% 
% % O = -[0 0.5];
% % landFunc = @(x,y) ((1-(x-O(1))).^2)+(100*(((y-O(2))-((x-O(1)).^2)).^2));        % Rosenbrock Function from CEC 2005
% % land_func_name = 'Rosenbrock';
% 
% % landFunc = @(x,y) sin(4*pi*x) + sin(4*pi*y);        % sine shape
% % land_func_name = 'sine';
% landFunc = @(x,y) sin(3*pi*x) + sin(3*pi*y);        % sine_2 shape
% land_func_name = 'sine2';
% 
% landFunc = @(x,y) sin(2*pi*x) + sin(2*pi*y);        % sine_2 shape
% land_func_name = 'sine2_2pix';

peakz = @(xx,yy) 3*(1-xx).^2.*exp(-(xx.^2) - (yy+1).^2) - 10*(xx/5 - xx.^3 - yy.^5).*exp(-xx.^2-yy.^2) - 1/3*exp(-(xx+1).^2 - yy.^2) ;
% landFunc = @(x,y) peakz(2*x,2*y) ;        % peak function
landFunc = @(x,y) peakz(x,y) ;        % peak function
land_func_name = 'peakz_x_y';

% sig = 0.1; mu = 0.2;
% gausFun = @(r) 1/sig/sqrt(2*pi)*exp(-0.5*((r-mu)/sig).^2);
% landFunc = @(x,y) gausFun(sqrt(x.^2+y.^2));
% land_func_name = 'donut';

% k = 2*pi;
% landFunc = @(x,y) 2*((x) + (y)).^1 + sin(k*x) + sin(k*y);%

ArenaScale = 1.0;
landBounds = ArenaScale*[-1 1; -1 1]; %
initBounds = landBounds;

n_tiles = 100;
x = linspace(landBounds(1,1),landBounds(1,2),n_tiles);
y = linspace(landBounds(2,1),landBounds(2,2),n_tiles);
[xx,yy] = meshgrid(x,y);
zz = landFunc(xx,yy);
% landFunc = @(x,y)landFunc(x,y);
landFunc = @(x,y) (landFunc(x,y)-min(zz(:)))/(max(zz(:))-min(zz(:)));
zz = landFunc(xx,yy);
zEnv = mean(zz(:));


% sigma = 0.01;
% z_noise = sigma*randn(size(zz));
% zz = zz + z_noise;

colormap gray
surf(xx,yy,zz,'EdgeAlpha',0,'FaceColor','flat')
view(0,90)
axis square
hold on
contour3(xx,yy,zz ,'Color',[238, 255, 150]/255,'LineWidth',1.25)
contour3(xx,yy,zz,[zEnv,zEnv],'color',[13, 182, 255]/255,'LineWidth',4)

set(gca,'Color','w')
set(gcf,'Color','w')
axis off

savefig(land_func_name+"_contours.fig")
exportgraphics(fig,land_func_name+"_contours.png",'resolution',300)
