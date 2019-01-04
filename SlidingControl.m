%%--------------------------
%----- SLIDING CONTROL------
%%--------------------------
clear;clc;
%% State space of system
a = 2;
b = 1;
c1 = 1;
c2 = 1;
Mu = 0.5;
t = (0:0.01:10)';
dmax = 0.9;
d = dmax*sin(628*t);

% initial state values
init_cond = [1 1]';

A = [a 0;1 0];
B = [b;0];
D = [1;0];
Cy = [1 0;0 1];

%% Sliding surface polynomial coefficients
surf_C = [c1;c2];

% surface stability check
if (-1)*surf_C(2)/surf_C(1) < 0
    disp('Stable');
else
    disp('Roots of surface polynomial is on closed right plane, unstable/marginally stable');
end

%% Plot of Phase Portrait and surface line with sliding control (signed function)
surfline = [(-1)*(-10:0.1:10)' (-10:0.1:10)'];
figure;
hold on;
plot(surfline(:,1),surfline(:,2),'-.r');
time_inc = 1e-3;

% pools of initial state values whose convergence towards zero to be
% investigated and plotted in phase portrait plane
x1_init = [-2 -1.5 -1 -0.5 0.5 1 1.5 2];
x2_init = x1_init;

for i=1:length(x1_init)
    x = [x1_init(i);x2_init(i)];
    xdot = x;
    surf_val = surf_C'*x;
    x_traj = x';
    t = 0;
    
    % looping to terminate the dynamics after reaching steady state
    % condition/reaching x-->0 close to zero, ie : [x1 x2] ~ [0 0]
    while sum(abs(x))>1e-3
        t = t + time_inc;
        u_con = (-1)*(surf_C'*A*x + Mu*sign(surf_val)+sign(surf_val)*surf_C'*D*dmax);
        u_con = u_con/(surf_C'*B);
        
        %state space and difference equation
        xdot = A*x + B*u_con + D*dmax*sin(628*t);
        x = x + xdot*time_inc;
        
        %evolution/movement of state values
        x_traj = [x_traj;x'];
        
        %movement of sliding polynomial values, ie: p1x1+p2x2
        surf_val = surf_C'*x;
        sumtrack = sum(abs(x)) % tracking convergence of state values towards zero/below acceptable small value (epsilon)
    end
    
    space = 200; %trajectory spacing to not dense the phase portrait arrows
    
    %getting direction of convergence
    p1 = x_traj(:,1);
    p2 = x_traj(:,2);
    pquiv1 = p1(1:space:end);
    pquiv2 = p2(1:space:end);
    dpquiv1 = gradient(pquiv1);
    dpquiv2 = gradient(pquiv2);
    
    
    plot(pquiv1,pquiv2,'-.b');
    quiver(pquiv1,pquiv2,dpquiv1,dpquiv2,'AutoScaleFactor',1.3,'MaxHeadSize', 1);
end
hold off;
grid;
axis([2*min(x1_init)-0.1 2*max(x1_init)+0.1 2*min(x1_init)-0.1 2*max(x1_init)+0.1]);
xlabel('X_1');
ylabel('X_2');
title('Phase portrait with various starting initial state values - signed sliding control');

%% Plot of Phase Portrait and surface line with smooth sliding control (saturated function)
surfline = [(-1)*(-10:0.1:10)' (-10:0.1:10)'];
figure;
hold on;
plot(surfline(:,1),surfline(:,2),'-.r');
time_inc = 1e-3;
epsilon = 0.01;

% pools of initial state values whose convergence towards zero to be
% investigated and plotted in phase portrait plane
x1_init = [-2 -1.5 -1 -0.5 0.5 1 1.5 2];
x2_init = x1_init;

for i=1:length(x1_init)
    x = [x1_init(i);x2_init(i)];
    xdot = x;
    surf_val = surf_C'*x;
    x_traj = x';
    t = 0;
    
    % looping to terminate the dynamics after reaching steady state
    % condition/reaching x-->0 close to zero, ie : [x1 x2] ~ [0 0]
    while sum(abs(x))>1e-3
        t = t + time_inc;
        u_con = (-1)*(surf_C'*A*x + Mu*sat(surf_C,x,epsilon)+sat(surf_C,x,epsilon)*surf_C'*D*dmax);
        u_con = u_con/(surf_C'*B);
        
        %state space and difference equation
        xdot = A*x + B*u_con + D*dmax*sin(628*t);
        x = x + xdot*time_inc;
        
        %evolution/movement of state values
        x_traj = [x_traj;x'];
        
        %movement of sliding polynomial values, ie: p1x1+p2x2
        surf_val = surf_C'*x;
        sumtrack = sum(abs(x)) % tracking convergence of state values towards zero/below acceptable small value (epsilon)
    end
    
    space = 200; %trajectory spacing to not dense the phase portrait arrows
    
    %getting direction of convergence
    p1 = x_traj(:,1);
    p2 = x_traj(:,2);
    pquiv1 = p1(1:space:end);
    pquiv2 = p2(1:space:end);
    dpquiv1 = gradient(pquiv1);
    dpquiv2 = gradient(pquiv2);
    
    plot(pquiv1,pquiv2,'-.b');
    quiver(pquiv1,pquiv2,dpquiv1,dpquiv2,'AutoScaleFactor',1.3,'MaxHeadSize', 1);
end
hold off;
grid;
axis([2*min(x1_init)-0.1 2*max(x1_init)+0.1 2*min(x1_init)-0.1 2*max(x1_init)+0.1]);
xlabel('X_1');
ylabel('X_2');
title('Phase portrait with various starting initial state values - saturated sliding control');

%% modularized function for surface saturation
function satval = sat(surfcoeff,state,epsilon)
    if surfcoeff'*state>epsilon
        satval = 1;
    elseif and((surfcoeff'*state>=(-1)*epsilon),(surfcoeff'*state<=epsilon))
        satval = surfcoeff'*state/epsilon;
    else
        satval = -1;
    end
end