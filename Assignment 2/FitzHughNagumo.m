close all
clear all

% Defining Variables
a = 0.5;
b = 0.1;
r = 0.1;

% Defining the FitzHugh-Nagumo equations
fv = @(v) v*(a-v)*(v-1);
dv_dt = @(v, w, I) fv(v) - w + I;
dw_dt = @(v, w) b*v - r*w;

v = -0.5 : 0.01 : 1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1 --------------------------------------------------
dt = 0.1;
tt = 100;
Iext = 0;
[vnc, wnc] = nullClines(v,Iext, b,r,fv,1);

v01 = 0.4;
v02 = 0.6;
w0 = 0;
t = 0.1:dt:tt;
[vhist1,whist1] = eulerInt(v01,w0,Iext,dt,tt, dv_dt,dw_dt);
[vhist2,whist2] = eulerInt(v02,w0,Iext,dt,tt, dv_dt,dw_dt);
plotfw(vhist1,whist1,t,0.7,2,1,Iext)
plotfw(vhist1,whist1,t,0.7,3,2,Iext)

v0 = -0.5:1:1.5;
w0 = -0.5:0.2:1.5;
for i=1:length(v0)
    for j=1:length(w0)
        [vhist3,whist3] = eulerInt(v0(i),w0(j),Iext,dt,tt, dv_dt,dw_dt);
        figure(1)
        plot(vhist3,whist3,'Color','#21DECC')
    end
end
xlim([-0.5 1.5]);

plot(vhist1,whist1,'Color','#6600FF',LineWidth=1.5)
plot(vhist2,whist2,'Color','#FF80ED',LineWidth=1.5)

% Case 2 --------------------------------------------------
dt = 0.1;
tt = 100;
Iext = 0.5;
[vnc, wnc] = nullClines(v,Iext, b,r,fv,4);

v0 = 0.4;
w0 = 0;
t = 0.1:dt:tt;
[vhist1,whist1] = eulerInt(v0,w0,Iext,dt,tt, dv_dt,dw_dt);
plotfw(vhist1,whist1,t,1.6,5,3,Iext)

v01 = 0.4;
w01 = 0.5;
[vhist2,whist2] = eulerInt(v01,w01,Iext,dt,tt, dv_dt,dw_dt);
figure(4)
plot(vhist2,whist2)

v0 = -0.5:1:1.5;
w0 = -0.5:0.2:1.5;
for i=1:length(v0)
    for j=1:length(w0)
        [vhist2,whist2] = eulerInt(v0(i),w0(j),Iext,dt,tt, dv_dt,dw_dt);
        figure(4)
        plot(vhist2,whist2,'Color','#21DECC')
    end
end
xlim([-0.5 1.5]);

% Case 3 --------------------------------------------------
dt = 0.1;
tt = 100;
Iext = 1;
[vnc, wnc] = nullClines(v,Iext, b,r,fv,6);

v0 = 0.4;
w0 = 0;
t = 0.1:dt:tt;
[vhist1,whist1] = eulerInt(v0,w0,Iext,dt,tt, dv_dt,dw_dt);
plotfw(vhist1,whist1,t,2,7,3,Iext)

v0 = -0.5:1:1.5;
w0 = -0.5:0.2:1.5;
for i=1:length(v0)
    for j=1:length(w0)
        [vhist2,whist2] = eulerInt(v0(i),w0(j),Iext,dt,tt, dv_dt,dw_dt);
        figure(6)
        plot(vhist2,whist2,'Color','#21DECC')
    end
end
xlim([-0.5 1.5]);

% Case 4 --------------------------------------------------
dt = 0.1;
tt = 100;
Iext = 0.05;
b = 0.01;
[vnc, wnc] = nullClines(v,Iext, b,r,fv,8);
ylim([-0.1 0.2]);

v0 = 0.4;
w0 = 0;
t = 0.1:dt:tt;
[vhist1,whist1] = eulerInt(v0,w0,Iext,dt,tt, dv_dt,dw_dt);
plotfw(vhist1,whist1,t,2,9,3,Iext)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for plotting nullclines and trajectories
function [vnc, wnc] = nullClines(v,Iext, b,r,fv,fig)
    iter = 1;
    for i = v
        vnc(iter) = fv(i) + Iext;
        wnc(iter) = b*i/r;
        iter = iter + 1;
    end
    figure(fig);
    plot(v,vnc,'b',v,wnc,'r',LineWidth=1.75)
    title(['Phase plot for I_{ext} = ', num2str(Iext)])
    xlabel('v','FontWeight','bold')
    ylabel('w','FontWeight','bold')
    text(0.5,Iext,'a','HorizontalAlignment','center','VerticalAlignment','top','Color','b')
    grid on
    hold on
    ylim([-0.5 1.5]);
end

function [vhist,whist] = eulerInt(v0,w0,Iext,dt,tt, dv_dt,dw_dt)
    niter = tt/dt;
    for i = 1:niter
        v0 = v0 + dv_dt(v0,w0,Iext)*dt; 
        w0 = w0 + dw_dt(v0,w0)*dt;
        vhist(i) = v0;
        whist(i) = w0;
    end
end

function [] = plotfw(vhist,whist,t,ymax,fig,subcase,Iext)
    s = {'I_{ext}=%s, Subcase 1: V(0) < a and W(0) = 0','I_{ext}=%s, Subcase 1: V(0) > a and W(0) = 0','I_{ext}=%s, V(0) < a and W(0) = 0'};
    figure(fig)
    plot(t,vhist,t,whist,LineWidth=1.75)
    ss = string(s(subcase));
    title(sprintf(ss,num2str(Iext)))
    legend('V(t)','W(t)')
    ylabel('V(t) and W(t)','FontWeight','bold')
    xlabel('t','FontWeight','bold')
    ylim([-0.2,ymax]);
    grid on
end