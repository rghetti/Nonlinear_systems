%% Settings

clear all
close all
clc
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex'); 

%% Data

data.alpha = 0.25;
data.beta = 0.5;
data.gamma_s = 0.35;
data.gamma_t = 0.25;
data.delta = 0.25;
data.epsilon = 0.65;
data.tau = 0.25;


%% Section 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that, in the following function, t is a scalar, while x is a 3D
% vector, whose coordinates can be accessed as follows: x(1), x(2), x(3).

u = @(t,x) x(1)-data.alpha/data.beta*x(2)+1/data.delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stop here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LVM = @(t,x) [x(1)*(1-x(1)-data.alpha*x(2)-data.gamma_s*x(3));
    x(2)*(1-x(2)-data.beta*x(1)-data.gamma_t*x(3));
    data.delta*x(3)*(-data.alpha/data.beta*x(2)+x(1))];

LVM_FL1 = @(t,x) [x(1)*(1-x(1)-data.alpha*x(2)-data.gamma_s*x(3));
    x(2)*(1-x(2)-data.beta*x(1)-data.gamma_t*x(3));
    data.delta*x(3)*(-data.alpha/data.beta*x(2)+x(1) - u(t,x))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [1,1,2];
t = linspace(0,200,200);

[t,y] = ode45(LVM,t,x0);
figure(1)
plot(t,y(:,1),'LineWidth',1.35)
hold on
plot(t,y(:,2),'LineWidth',1.35)
plot(t,y(:,3),'LineWidth',1.35)
xlabel("$t$","FontSize",35,'Interpreter','latex')
title("Simulation of the open-loop system",'FontSize',30,'Interpreter','latex')
legend("$S(t)$","$T(t)$","$A(t)$","Fontsize",30,'Interpreter','latex','location','east')
set(figure(1),'Position', get(figure(1), 'Position')+[-320, -200, 700, 400]); 
grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = linspace(0,20,200);
[t,y] = ode45(LVM_FL1,t,x0);
figure(2)
plot(t,y(:,1),'LineWidth',1.35)
hold on
plot(t,y(:,2),'LineWidth',1.35)
plot(t,y(:,3),'LineWidth',1.35)
xlabel("$t$","FontSize",35,'Interpreter','latex')
title("Feedback Linearisation Controller $u$ - Simulation of the System",'FontSize',30,'Interpreter','latex')
legend("$S(t)$","$T(t)$","$A(t)$","Fontsize",30,'Interpreter','latex','location','northeast')
set(figure(2),'Position', get(figure(2), 'Position')+[-320, -200, 700, 400]); 
grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = arrayfun(@(i) u(1,y(i, :)'), 1:size(y, 1))';
figure(3)
plot(t,U,'LineWidth',1.35)
hold on
xlabel("$t$","FontSize",35,'Interpreter','latex')
title("Feedback Linearisation Controller $u$ - Control Effort",'FontSize',30,'Interpreter','latex')
legend("$U(t)$",'Fontsize',30,'Interpreter','latex','location','east')
set(figure(3),'Position', get(figure(3), 'Position')+[-320, -200, 700, 400]); 
grid on
hold off



%% Section 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = @(t,x) 1/(x(3)^2*x(4)*data.delta*data.epsilon)*(-data.delta*data.epsilon*x(4)*x(3)*(1-x(4)+data.tau*x(3)) ...
            +data.delta^2*x(3)*(-data.alpha/data.beta*x(2)+x(1)-x(4))^2-data.alpha*data.delta/data.beta*x(3)*x(4) ...
            *(1-x(2)-data.beta*x(1)-data.gamma_t*x(3))+data.delta*x(3)*x(1)*(1-x(1)-data.alpha*x(2)-data.gamma_s*x(3))+x(3) ...
            +2*(data.delta*x(3)*(-data.alpha/data.beta*x(2)+x(1)-x(4))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stop here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LVM_FL2 = @(t,x) [x(1)*(1-x(1)-data.alpha*x(2)-data.gamma_s*x(3));
    x(2)*(1-x(2)-data.beta*x(1)-data.gamma_t*x(3));
    data.delta*x(3)*(-data.alpha/data.beta*x(2)+x(1)-x(4));
    data.epsilon*x(4)*(1-x(4)+data.tau*x(3) + x(3)*w(t,x) )];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [1,1,1,0.01]';
t = linspace(0,40,10^2);
[t,y] = ode45(LVM_FL2,t,x0);
figure(2)
plot(t,y(:,1),'LineWidth',1.35)
hold on
plot(t,y(:,2),'LineWidth',1.35)
plot(t,y(:,3),'LineWidth',1.35)
plot(t,y(:,4),'LineWidth',1.35)
xlabel("$t$","FontSize",35,'Interpreter','latex')
title("Feedback Linearisation Controller $w$ - Simulation of the System",'FontSize',30,'Interpreter','latex')
legend("$S(t)$","$T(t)$","$A(t)$","$P(t)$","Fontsize",30,'Interpreter','latex','location','northeast')
set(figure(2),'Position', get(figure(2), 'Position')+[-320, -200, 700, 400]); 
grid on
hold off
print -depsc 'k.eps'; 




