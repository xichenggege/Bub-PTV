%##########################################################################
% Bubble tracking post processing
% compute 3D trajetory and velocity of the data from image-processing
%##########################################################################
% Demo version v1.0
% 2022-06-12
% © Xicheng Wang & Yun Feng
%==========================================================================
clc
clear all
close all

% %% Post processing global  setup
run('InitialSetup_SEFW7_V1');   % run setup file
addpath('.\Functions');

figureFolder = strcat('Figures\post-processing\');
dataFolder   = strcat('Data\post-processing\');
mkdir(strcat('Figures\post-processing\'));

Vlimit = 1.2;

% Note this info was derived from calibration test, i.e. with checkboard
% Offset of global coordinate system [mm] 
% Thus, currently, the coordinate should be revised as below
% to move origin to orifice local coordinate
L_off = (32.6718+32.4422)*0.5-20;  % 10mm is checkboard size, 1 block length   
H_off = -(10.3917+20.3918)*0.5;  
D_off = (661.4654+661.2636)*0.5;

% grid info
Lmin = 20;   Lmax = 320;
Dmin = -40;  Dmax = 40;
Hmin = -80;  Hmax = 80;

nL = 61; nD = 17; nH = 33;
L  = linspace (Lmin,Lmax,nL);
D  = linspace(Dmin,Dmax,nD);
H  = linspace(Hmin,Hmax,nH);

%% Load bubble information
% read bubble information, for tutorial purpose
file  = load(strcat(dataFolder,'cbbls'));
cbbls = file.cbbls;
file  = load(strcat(dataFolder,'Traj'));
Traj  = file.Traj;

% Trajectory offset 
for kk=1:numel(Traj)
    Traj{kk}.pos(:,1) = Traj{kk}.pos(:,1) - L_off;
    Traj{kk}.pos(:,2) = Traj{kk}.pos(:,2) - H_off;
    Traj{kk}.pos(:,3) = Traj{kk}.pos(:,3) - D_off;
end
griddata_available = false;

%% Frame length 
frameLength = [];
for kk=1:numel(cbbls)
    frameLength(kk) = cbbls{kk}.LastFrame-cbbls{kk}.FirstFrame;
end
set(0,'DefaultFigureVisible', 'on');
FigHandle = figure; hold on
set(FigHandle, 'Position', [50, 300, 600, 510]);

h = histogram(frameLength);
set(gca,'yscale','log');

ylabel('Frequency','interpreter','latex','fontsize',16,'fontname','times');
xlabel('Bubble frame length','interpreter','latex','fontsize',16,'fontname','times');
set(gca,'fontsize',16,'fontname','times');
set(gcf,'PaperPositionMode','auto');
saveas(FigHandle,strcat(figureFolder,'frameLength.png'));
    
%% plot 3D trajectory
file = load('rainbow');  % load colors
rainbow = file.rainbow;  

set(0,'DefaultFigureVisible', 'on');
FigHandle = figure;
set(FigHandle, 'Position', [150, 200, 800, 600]);
plotFreq = 2;
kk = 1;

for j = 1:round(numel(Traj)/plotFreq)
    % collect 'plotFreq' to draw together 
    coord = [];
    for i=1:plotFreq
        if kk<=numel(Traj)
            coord = [coord; Traj{kk}.pos(:,:)];
            kk=kk+1;
        end
    end
    % display color 
    if j<=length(rainbow)
        co = rainbow(j,:);
    else
        co = rand(1,3);
    end
    % plot
    plot3(coord(:,1),coord(:,3),coord(:,2),'.','color',co); hold on
end

set(gca,'ztick',[-80:40:80]);
set(gca,'ytick',[-80:40:80]);
set(gca,'xtick',[50:50:450]);

set(gca,'fontsize',16,'fontname','times');
zlabel('$z\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
ylab=ylabel('$y\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
xlab=xlabel('$x\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
%     set(xlab,'Rotation',23);
set(ylab,'Rotation',-30);

axis equal
view([-29.483987269116692,23.908035932214432]);
grid on

set(gca,'zlim',[-80 80]);
set(gca,'xlim',[20 320]);
set(gca,'ylim',[-42 42]);

exportgraphics(FigHandle,strcat(figureFolder,'trajectories.png'),'Resolution',450);


%% Velocity calculation
step         = 10;               % velocity calculation per 'step' intervals (time lag, better to be odd number)
FPS          = 6300;

% Velocity calculation
Traj         = vel_calculate(Traj,step,FPS);

% filter out too short trajectory
for kk=1:numel(Traj)
    Traj{kk}.LastFrame = NumFrames;
    flength(kk) = Traj{kk}.LastFrame - Traj{kk}.FirstFrame;
end
ind_empty = find(flength<=5);
Traj(ind_empty) = [];

%% Distribution of PTV results on grid
fitMethod  = char('natural');   % data interpolation method
[xgrid, zgrid, ygrid] = meshgrid(L,H,D); % 3D grid

world_time   = [1:1:NumFrames]*1/FPS;
dt           = 1/FPS;

if ~griddata_available
%     3D grid 
    Ffield3d   = btv2grid3(Traj,world_time,xgrid,L,H,D,FPS,fitMethod);   % Frame Field data
%     Save grid data
    save(strcat(dataFolder,'grid3'),'Ffield3d');
else
    file = load(strcat(dataFolder,'grid3'));
    Ffield3d = file.Ffield3d;
end

%-----------------------------------
% 3D grid data
%-----------------------------------
[Lgrid, Hgrid, Dgrid] = meshgrid(L,H,D); % 2D grid
bbcount3 = 0*Lgrid;
for f=1:NumFrames
    bbcount3 = bbcount3 + Ffield3d{f}.nbb;
end

l = Lgrid(:);
h = Hgrid(:);
d = Dgrid(:);
nbb = bbcount3(:);

% filter out too less data point
nbbfilt = 5;
nbb (nbb<nbbfilt) = nan;

% PLOT
set(0,'DefaultFigureVisible', 'on');
FigHandle = figure;
set(FigHandle, 'Position', [50, 300, 720, 500]);
scatter3(l,d,h,10,nbb,'filled');

colormap jet;
cb=colorbar('position',[0.92, 0.25, 0.025,0.6]);
set(cb,'fontsize',14,'fontname','times');                   
caxis([20 200]);

set(gca,'fontsize',18,'fontname','times');

set(gca,'ztick',[-80:40:80]);
set(gca,'ytick',[-80:20:80]);
set(gca,'xtick',[50:100:450]);

zlabel('$z\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
ylab=ylabel('$y\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
xlab=xlabel('$x\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
set(ylab,'Rotation',-10);

% axis equal
view([-54.615935942507683,18.012281297363188]);
grid on

set(gca,'zlim',[-80 80]);
set(gca,'xlim',[20 320]);
set(gca,'ylim',[-42 42]);

set(FigHandle ,'PaperPositionMode','auto'); 
exportgraphics(FigHandle,strcat(figureFolder,'BubbleDensity3d.png'),'Resolution',450);

%% Result averaging, PIV similar post-processing
for f=1:numel(Ffield3d)
    UField{f}    = Ffield3d{f}.U./ Ffield3d{f}.nbb;
end
% Filter by 1.5 std
stdevMultiplier = 1.5; % remove data out of 1.5 sigma range;
[U_median, U_avg, ~,~] = fieldAveraging(UField,'stdFilter',stdevMultiplier);

%% Flow reconstruction by surface fitting and scatter interpolation
% roughly estimation, should be confirmed
l  = linspace (Lmin,Lmax,nL);
d  = linspace(Dmin,Dmax,nD);
h  = linspace(Hmin,Hmax,nH);

[xgrid, ygrid, zgrid] = meshgrid(l,d,h); % 3D grid

xslice = [0.15 0.25]*1e3;
yslice = [0.00]*1e3;
zslice = [-0.01]*1e3;

%% Raw averaged 3D data
ydata = permute(U_avg,[3 2 1])/1e3;
FigHandle = figure;
set(FigHandle, 'Position', [50, 300, 800, 600]);
slice(xgrid,ygrid,zgrid,ydata,xslice,yslice,zslice); 

colormap jet;
caxis([0 Vlimit]);    % [m/s]
cb=colorbar('position',[0.92, 0.33, 0.025,0.55]);
set(cb,'YTick',0:0.3:1.2,'fontsize',14,'fontname','times');  
set(cb,'yTickLabel',num2str(get(cb,'yTick')','%.1f'))

set(gca,'ztick',[-80:40:80]);
set(gca,'ytick',[-80:40:80]);
set(gca,'xtick',[50:50:450]);

set(gca,'fontsize',16,'fontname','times');

zlabel('$z\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
ylab=ylabel('$y\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
xlab=xlabel('$x\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
set(ylab,'Rotation',-30);

axis equal
view([-29.483987269116692,23.908035932214432]);
grid on

set(gca,'zlim',[-80 80]);
set(gca,'xlim',[20 320]);
set(gca,'ylim',[-42 42]);

text (340,0,95,['$U\ [m/s]$'],'interpreter','latex','fontsize',16);
exportgraphics(FigHandle,strcat(figureFolder,'Uavg3D_raw.png'),'Resolution',450);


% Flow reconstruction by averaged 3D data
myfun = @(paras,x)(paras(1)./(x(:,1)-paras(3))).*exp(-paras(2)*(x(:,2).^2+x(:,3).^2)./((x(:,1)-paras(3)).^2));
% from single-phase round jet,  U0 B d and K
paras0 = [0.14 77 -0.025];
up = [0.5 150 0];
lb = [0.01 5 -0.1];

%% fitting surface 
% Convert to 1 column
ydata = permute(U_avg,[3 2 1])/1e3;  % mm/s to m/s
yfit  = ydata(:);
xfit  = [xgrid(:) ygrid(:) zgrid(:)]/1e3;
xrec  = xfit;

% remove nan
idx = isnan(yfit);
xfit(idx,:) = [];
yfit(idx) = [];

% remove negative
idx = find(yfit<0);
xfit(idx,:) = [];
yfit(idx) = [];

% Fitting buoyant jet equation
if numel(yfit)>10
    [paras,resnorm,residual] = lsqcurvefit(myfun,paras0,xfit,yfit,lb,up);
    x = xrec;
    ypredict = (paras(1)./(x(:,1)-paras(3))).*exp(-paras(2)*(x(:,2).^2+x(:,3).^2)./((x(:,1)-paras(3)).^2));
else
    paras = nan*ones(1,numel(paras0));
    x = xrec;
    ypredict = ones(numel(xrec(:,1)),1)*nan;
end

% Reconstruction
Urecon = reshape(ypredict,size(xgrid));
% plot - Reconstruction based on raw averaged 3D data
FigHandle = figure;
set(FigHandle, 'Position', [50, 300, 800, 600]);
slice(xgrid,ygrid,zgrid,Urecon,xslice,yslice,zslice); 
shading interp;

colormap jet;
caxis([0 Vlimit]);    % [m/s]
cb=colorbar('position',[0.92, 0.33, 0.025,0.55]);
set(cb,'YTick',0:0.3:1.2,'fontsize',14,'fontname','times');  
set(cb,'yTickLabel',num2str(get(cb,'yTick')','%.1f'))

set(gca,'ztick',[-80:40:80]);
set(gca,'ytick',[-80:40:80]);
set(gca,'xtick',[50:50:450]);

set(gca,'fontsize',16,'fontname','times');

zlabel('$z\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
ylab=ylabel('$y\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
xlab=xlabel('$x\ [mm]$','interpreter','latex','fontsize',20,'fontname','times');
set(ylab,'Rotation',-30);

axis equal
view([-29.483987269116692,23.908035932214432]);
grid on

set(gca,'zlim',[-80 80]);
set(gca,'xlim',[20 320]);
set(gca,'ylim',[-42 42]);

text (340,0,95,['$U\ [m/s]$'],'interpreter','latex','fontsize',16);
exportgraphics(FigHandle,strcat(figureFolder,'Uavg3D_reconstructed.png'),'Resolution',450);

%% 2D axial compressed -> assume axis-symmetry
% 2D coordinates
L  = linspace (Lmin,Lmax,nL);
R  = linspace(Hmin,Hmax,nH);
[lgrid, rgrid] = meshgrid(L,R);
Traj_off = Traj;
for kk=1:numel(Traj_off)
    Traj_off{kk}.Cellcenter(:,1) = Traj_off{kk}.Cellcenter(:,1);
    Traj_off{kk}.Cellcenter(:,2) = Traj_off{kk}.Cellcenter(:,2);
    Traj_off{kk}.Cellcenter(:,3) = Traj_off{kk}.Cellcenter(:,3);
end

% 2D grid 
Ffield2d = btv2grid2(Traj_off,world_time,lgrid,L,R,FPS);   % Frame Field data
save(strcat(dataFolder,'grid2'),'Ffield2d');

% velocity averaging 2D
for f=1:numel(Ffield2d)
    UField2{f}    = Ffield2d{f}.U_avg2./ Ffield2d{f}.nbb2;
end
% Filter by 3.0 std
stdevMultiplier = 3; % remove data out of 1.5 sigma range;
[~, U_avg2, ~,~]  = fieldAveraging(UField2,'stdFilter',stdevMultiplier);

U_avg2    = U_avg2./1e3;
% 2D Plot and compare with reconstructed and CFD
FigHandle = figure;
set(FigHandle, 'Position', [50, 300, 600, 400]);
pcolor(lgrid,rgrid,U_avg2);
colormap jet;
shading interp;
cb=colorbar;
set(cb,'YTick',0:0.3:1.2,'fontsize',14,'fontname','times');  
set(cb,'yTickLabel',num2str(get(cb,'yTick')','%.1f'));
caxis([0 Vlimit]);    % [m/s]

set(gca,'ytick',[-120:40:120]);
set(gca,'xtick',[0:50:300]);

text (320,115,['$U\ [m/s]$'],'interpreter','latex','fontsize',14);     
% Plot sparger location
xlabel('$x\ [mm]$ ','interpreter','latex');
ylabel('$r\ [mm]$ ','interpreter','latex');
set(gca,'fontsize',18,'fontname','times'); 
set(FigHandle ,'PaperPositionMode','auto'); grid off;
axis equal

set(gca,'ylim',[-100 100]);
set(gca,'xlim',[20 300]);

exportgraphics(FigHandle,strcat(figureFolder,'Uavg_axis_sym_V1.png'),'Resolution',450);


%% 2D Plot Reconstructed
% Convert to 1 column
yfit  = U_avg2(:);
xfit  = [lgrid(:) rgrid(:)]/1e3;
xrec  = xfit;

% remove nan
idx = isnan(yfit);
xfit(idx,:) = [];
yfit(idx) = [];

% remove negative
idx = find(yfit<0);
xfit(idx,:) = [];
yfit(idx) = [];

myfun = @(paras,x)(paras(1)./(x(:,1)-paras(3))).*exp(-paras(2)*(x(:,2).^2)./((x(:,1)-paras(3)).^2));
% from single-phase round jet,  U0 B d and K
paras0 = [0.14 77 -0.025];
up = [0.5 100 0];
lb = [0.01 5 -0.1];

% Fitting buoyant jet equation
[paras2d,resnorm,residual] = lsqcurvefit(myfun,paras0,xfit,yfit,lb,up);
x = xrec;
ypredict = (paras2d(1)./(x(:,1)-paras2d(3))).*exp(-paras2d(2)*(x(:,2).^2)./((x(:,1)-paras2d(3)).^2));

% Reconstruction
Urecon2 = reshape(ypredict,size(lgrid));

FigHandle = figure;
set(FigHandle, 'Position', [50, 300, 600, 400]);
pcolor(lgrid,rgrid,Urecon2);
colorbar;
colormap jet;
shading interp;
cb=colorbar;
set(cb,'YTick',0:0.3:1.2,'fontsize',14,'fontname','times');  
set(cb,'yTickLabel',num2str(get(cb,'yTick')','%.1f'));

caxis([0 Vlimit]);    % [m/s]

set(gca,'ytick',[-120:40:120]);
set(gca,'xtick',[0:50:300]);

text (320,115,['$U\ [m/s]$'],'interpreter','latex','fontsize',14);     
% Plot sparger location
xlabel('$x\ [mm]$ ','interpreter','latex');
ylabel('$r\ [mm]$ ','interpreter','latex');
set(gca,'fontsize',18,'fontname','times'); 
set(FigHandle ,'PaperPositionMode','auto'); grid off;
axis equal

set(gca,'ylim',[-100 100]);
set(gca,'xlim',[20 300]);

exportgraphics(FigHandle,strcat(figureFolder,'Urecon_axis_sym_V1.png'),'Resolution',450);
