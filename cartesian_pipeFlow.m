% Vorticity-Stream Function Formulation in Cartesian Coordinates
% Set up for pipe flow boundary conditions
% Ayaaz Yasin - Nov 02, 2024
clear; clc; close all;

%%% Problem Setup
rho     = 1;            %[kg/m^3],      mass density
mu      = 1;            %[Pa-s],        dynamic viscosity
uMax    = 1;            %[m/s],         inlet velocity
xMax    = 10;           %[m],           domain size in x
yMax    = 1;            %[m],           domain size in y

L       = yMax;         %[m],           characteristic length
Re      = rho*uMax*L/mu;%[-],           reynolds number
monitor = 0;

%%% Numerical Setup
nx      = 501;          %[-],           x-grid size
ny      = 51;           %[-],           y-grid size
minIter = 1e4;          %[-],           min iterations required
maxIter = 1e10;         %[-],           max iterations allowed
relax   = 0.1;          %[-],           vorticity relaxation factor
tol     = 1e-6;        %[-],           convergence tolerance
x = linspace(0,xMax,nx);    dx = x(2)-x(1);
y = linspace(0,yMax,ny);    dy = y(2)-y(1);

%%% Boundary Conditions (set for pipe flow. left inlet and right zero-gradient boundary.)
u_top   = 0;        u_left  = uMax;    u_right     = 0;    u_bottom    = 0;
v_top   = 0;        v_left  = 0;        v_right     = 0;    v_bottom    = 0;
s_top   = uMax*max(y);   s_left  = uMax*y;    s_right     = 0;    s_bottom    = min(y);


%%% Intialization
initial = 0;
%psi function and new
s = initial*ones(length(x), length(y));    
s(1,:) = s_left;   s(end,:) = s_right; s(:,1) = s_bottom;   s(:,end) = s_top;
sNew = s;

%u-velocity
u = initial*ones(length(x), length(y));
u(1,:) = u_left;   u(end,:) = u_right; u(:,1) = u_bottom;   u(:,end) = u_top;

%v-velocity
v = initial*ones(length(x), length(y));
v(1,:) = v_left;   v(end,:) = v_right; v(:,1) = v_bottom;   v(:,end) = v_top;

w = initial*ones(length(x), length(y));    wNew = w;   %omega and new
dpdx = nan(length(x), length(y));



%% Solution Loop
iter = 0;
while true
    iter = iter+1;
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            %Update psi
            term1 = (dy^2)*(s(i+1,j) + s(i-1,j));
            term2 = (dx^2)*(s(i,j+1) + s(i,j-1));
            term3 = w(i,j)*(dx^2)*(dy^2);
            denom = 2*(dx^2 + dy^2);
            check = (term1+term2+term3)/denom;
            sNew(i,j) = check;

            %Update omega
            term1 = (dy^2)*(w(i+1,j) + w(i-1,j));
            term2 = (dx^2)*(w(i,j+1) + w(i,j-1));
            A1 = (s(i,j+1)-s(i,j-1))/(2*dy);
            A2 = (w(i+1,j)-w(i-1,j))/(2*dx);
            A3 = (s(i+1,j)-s(i-1,j))/(2*dx);
            A4 = (w(i,j+1)-w(i,j-1))/(2*dy);
            A = A1*A2 - A3*A4;
            term3 = -Re*A*(dx^2)*(dy^2);
            denom = 2*(dx^2 + dy^2);
            check = (term1+term2+term3)/denom;
            wNew(i,j) = check;

            %Update velocities
            u(i,j) =  (s(i,j+1) - s(i,j-1))/(2*dy);
            v(i,j) = -(s(i+1,j) - s(i-1,j))/(2*dx);
        end
    end

    %Enforce Boundary Conditions
    wNew(:,end) = -(u(:,end) - u(:,end-1))/dy;      %top
    wNew(:,1)   = -(u(:,2) - u(:,1))/dy;            %bottom
    wNew(1,:)   =  (v(2,:) - v(1,:))/dx;            %left
    wNew(end,:) =  wNew(end-1,:);%(v(end,:) - v(end-1,:))/dx;      %right
    
    sNew(end,:) =  sNew(end-1,:);   %outlet zero-gradient boundary
    u(end,:) = u(end-1,:);
    v(end,:) = v(end-1,:);

    %check convergence
    pErr = abs(max(max(sNew - s)));
    wErr = abs(max(max(wNew - w)));
    if sum(sum(isnan(sNew)))>1 || sum(sum(isnan(wNew)))>1; fprintf('vorticit-streamfunction diverged\n'); break; end
    if pErr < tol && wErr < tol && iter > minIter; fprintf('vorticit-streamfunction converged\n'); break; end
    if iter > maxIter; fprintf('iteration maxed\n'); break; end

    if monitor == 1; figure(1); 
        semilogy(iter, pErr, 'b.', iter, wErr, 'r.'); hold on; legend('psi','omega','location','northeast')
        xlabel('Iteration'); ylabel('Error'); set(gca, 'FontSize',14, 'FontWeight', 'bold'); axis padded
    end
    fprintf('%0.4e\t\t%0.4e\n', pErr, wErr)

    %shift variables
    s = sNew;
    w = relax*wNew + (1-relax)*w;
end


%%% pressure gradient
for i = 2:length(x)-1
    for j = 2:length(y)-1
        term1 = (1/Re)*((u(i-1,j) - 2*u(i,j) + u(i+1,j))/(dx^2) + (u(i,j-1) - 2*u(i,j) + u(i,j+1))/(dy^2));
        term2 = - u(i,j)*(u(i+1,j) - u(i-1,j))/(2*dx);
        term3 = - v(i,j)*(u(i,j+1) - u(i,j+1))/(2*dy);
        dpdx(i,j) = term1 + term2 + term3;
    end
end

%%% pressure field
p = nan(length(x), length(y)); p(end, :) = 0;
for i = length(x)-1:-1:1
    p(i,:) = p(i+1,:) - dx*dpdx(i,:);
end

%% Post-processing
close all; 
[X,Y] = meshgrid(x,y);

fig1=figure; fig1.Position=[-2065,65,1449,920];
subplot(3,1,1); hold on
uMag = sqrt(u.^2 + v.^2);
contourf(X,Y,uMag', 70, 'EdgeColor','none'); colormap(jet), colorbar(gca, "southoutside");
axis equal; xlim([0,x(end)]);ylim([0,y(end)])
hold off; set(gca, 'FontSize',14, 'FontWeight', 'bold')
title(sprintf('Re = %i, velocity contour',Re))

subplot(3,1,2); hold on
nPts = 30; 
xstart = min(x) + (max(x)-min(x))*zeros(nPts,1);
ystart = linspace(min(y), max(y), nPts); %min(y) + (max(y)-min(y))*rand(nPts,1);
streamline(X,Y,u',v',xstart, ystart, [0.5, 10000])
axis equal; xlim([0,x(end)]);ylim([0,y(end)])
hold off; set(gca, 'FontSize',14, 'FontWeight', 'bold')
title('streamlines')

subplot(3,1,3); hold on
contourf(X,Y,p', 70, 'EdgeColor','none'); colormap(jet), colorbar(gca, "southoutside");
axis equal; xlim([0,x(end)]);ylim([0,y(end)])
hold off; set(gca, 'FontSize',14, 'FontWeight', 'bold')
title('pressure')

%%% velocity cross-section
fig2=figure; fig2.Position=[-585   563   560   420];
xpos = max(x); idx = find(x == xpos);
uMag1 = uMag';
plot(uMag1(:,idx), y, 'r-', 'LineWidth',2)
xlabel('velocity magnitude [m/s]'); ylabel('y-coordinate [m]')
axis padded; set(gca, 'FontSize',14, 'FontWeight', 'bold')