function plot_traj(n)

load("new_data.mat")
load("WIcover.mat")
load("IDindx.mat")

% Plot the whole trajectory of the nth deer
xll=-10124202; dx=30;  % lower left corner and pixel size in meters
yll=5190957;
nrow=6191;ncol= 5860;
x=[0:ncol-1]*dx;
y=[0:nrow-1]*dx;
[X,Y]=meshgrid(x,y); 

xs=new_data(IDindx(n,1):IDindx(n,2), 1)-xll; 
ys=new_data(IDindx(n,1):IDindx(n,2),2)-yll;
maxx=max(xs); minx=min(xs); lx=maxx-minx;  
maxy=max(ys); miny=min(ys); ly=maxy-miny;

nbuf=2; 
jcolmin=floor( (minx)/dx )+1-nbuf; 
irowmin=floor( (miny)/dx )+1-nbuf;
jcolmax=floor( (maxx)/dx )+1+nbuf;
irowmax=floor( (maxy)/dx )+1+nbuf;
indxx=[jcolmin:jcolmax]; indxy=[irowmin:irowmax];
cover_small=WIcover(indxy,indxx);
X_small=X(indxy,indxx); Y_small=Y(indxy,indxx); 

pcolor(X_small/1000,Y_small/1000,cover_small), shading flat , axis image
hold on, plot(xs/1000,ys/1000,'k'), hold off
xlabel('x (km)'), ylabel('y (km)')
title(['Deer ', num2str(n)])