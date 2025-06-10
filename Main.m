%task2
w1 = 1;
w2 = 0.05;
w3 = 100;
N=10;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
x=[0.9492 0.43]';
Xsim=zeros(2,Sim_time);
Xsim(:,1)=x;
Usim=zeros(1,Sim_time-1);
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim(:,i),N,xr1,w1,w2,w3,1);
        Usim(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim(i)=u;
        Xsim(:,i+1)=x;
    end
end
c=linspace(1,Sim_time,Sim_time)/5;
figure;
subplot(2, 1, 1);
plot(c,Xsim(2,:),'r','DisplayName','x_2')
hold on;
plot(c,Xsim(1,:),'b','DisplayName','x_1')
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'m','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'m','LineStyle','--','DisplayName','x_r2')
xlabel('Time (s)');
ylabel('Value');
title('System evolution')
legend;
grid on;

subplot(2,1,2);
plot(c(1:end-1),Usim,'b','DisplayName','Input')
legend;
xlabel('Time (s)');
ylabel('u(k)');
title('Optimal input')
grid on;

%% task 3

w1 = 1;
w2 = 0.05;
w3 = 100;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
%for N=5
N=5;
x=[0.9492 0.43]';
Xsim5=zeros(2,Sim_time);
Xsim5(:,1)=x;
Usim5=zeros(1,Sim_time-1);

for i=1:Sim_time-1
    disp(i)
    disp(N)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim5(:,i),N,xr1,w1,w2,w3,1);
        Usim5(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim5(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim5(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim5(i)=u;
        Xsim5(:,i+1)=x;
    end
end
%for N=10
N=10;
x=[0.9492 0.43]';
Xsim10=zeros(2,Sim_time);
Xsim10(:,1)=x;
Usim10=zeros(1,Sim_time-1);
for i=1:Sim_time-1
    disp(i)
    disp(N)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim10(:,i),N,xr1,w1,w2,w3,1);
        Usim10(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim10(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim10(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim10(i)=u;
        Xsim10(:,i+1)=x;
    end
end
%for N=20
N=20;
x=[0.9492 0.43]';
Xsim20=zeros(2,Sim_time);
Xsim20(:,1)=x;
Usim20=zeros(1,Sim_time-1);

for i=1:Sim_time-1
    disp(i)
    disp(N)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim20(:,i),N,xr1,w1,w2,w3,1);
        Usim20(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim20(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim20(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim20(i)=u;
        Xsim20(:,i+1)=x;
    end
end
c=linspace(1,Sim_time,Sim_time)'/5;

figure;
subplot(2,1,1);
plot(c,Xsim5(2,:),'r','DisplayName','Output5')
hold on
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'g','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'g','LineStyle','--','DisplayName','x_r2')
plot(c,Xsim10(2,:),'b','DisplayName','Output10')
plot(c,Xsim20(2,:),'k','DisplayName','Output20')
xlabel('Time(s)');
ylabel('x_2(k)');
legend;
title("Concentration level for different Ns")

subplot(2,1,2);
plot(c(1:end-1),Usim5,'r','DisplayName','Input5')
hold on
plot(c(1:end-1),Usim10,'b','DisplayName','Input10')
plot(c(1:end-1),Usim20,'k','DisplayName','Input20')
legend;
xlabel('Time(s)');
ylabel('u(k)');
title('Input intensity for different Ns')

%%
%task4
w1 = 1;
w2 = 0.05;
w3 = 10;
N=10;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
x=[0.9492 0.43]';
Xsim4=zeros(2,Sim_time);
Xsim4(:,1)=x;
Usim4=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    disp(w3)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim4(:,i),N,xr1,w1,w2,w3,1);
        Usim4(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim4(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim4(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim4(i)=u;
        Xsim4(:,i+1)=x;
    end
end

time10=toc;
%%the one from task2
w3=100;
x=[0.9492 0.43]';
Xsim2=zeros(2,Sim_time);
Xsim2(:,1)=x;
Usim2=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    disp(w3)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim2(:,i),N,xr1,w1,w2,w3,1);
        Usim2(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim2(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim2(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim2(i)=u;
        Xsim2(:,i+1)=x;
    end
end
time100=toc;
%plotting the two simulations in three differents plots
c=linspace(1,Sim_time,Sim_time)/5;
figure;
subplot(3, 1, 1);
plot(c,Xsim2(2,:),'r','DisplayName','100')
hold on;
plot(c,Xsim4(2,:),'b','DisplayName','10')
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'m','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'m','LineStyle','--','DisplayName','x_r2')
legend;
xlabel('Time(s)');
ylabel('Concentration level');
grid on

subplot(3,1,2)
plot(c,Xsim2(1,:),'r','DisplayName','100')
hold on
plot(c,Xsim4(1,:),'b','DisplayName','10')
legend;
xlabel('Time(s)');
ylabel('Temperature level');
grid on

subplot(3,1,3)
plot(c(1:end-1),Usim2,'r','DisplayName','100')
hold on
plot(c(1:end-1),Usim4,'b','DisplayName','10')
xlabel('Time (s)');
ylabel('u(k)');
legend;
grid on;
time100
time10

%%
%task5

w1 = 1;
w2 = 0.05;
w3 = 100;
N=10;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
x=[0.9492 0.43]';
Xsim2=zeros(2,Sim_time);
Xsim2(:,1)=x;
Usim2=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim2(:,i),N,xr1,w1,w2,w3,1);
        Usim2(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim2(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim2(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim2(i)=u;
        Xsim2(:,i+1)=x;
    end
end
%changing the objective function
x=[0.9492 0.43]';
Xsim5=zeros(2,Sim_time);
Xsim5(:,1)=x;
Usim5=zeros(1,Sim_time-1);
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim5(:,i),N,xr1,w1,w2,w3,2);
        Usim5(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim5(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim5(:,i),N,xr2,w1,w2,w3,2);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim5(i)=u;
        Xsim5(:,i+1)=x;
    end
end

c=linspace(1,Sim_time,Sim_time)/5;
figure;
subplot(3, 1, 1);
plot(c,Xsim2(2,:),'r','DisplayName','Task2')
hold on;
plot(c,Xsim5(2,:),'b','DisplayName','Task5')
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'m','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'m','LineStyle','--','DisplayName','x_r2')
legend;
xlabel('Time(s)');
ylabel('Concentration level');
grid on

subplot(3,1,2)
plot(c,Xsim2(1,:),'r','DisplayName','Task2')
hold on
plot(c,Xsim5(1,:),'b','DisplayName','Task5')
legend;
xlabel('Time(s)');
ylabel('Temperature level');
grid on

subplot(3,1,3)
plot(c(1:end-1),Usim2,'r','DisplayName','Task2')
hold on
plot(c(1:end-1),Usim5,'b','DisplayName','Task5')
xlabel('Time (s)');
ylabel('u(k)');
legend;
grid on;
%%
%task 6: same simulation but adding a noise to the input of NMPC and
%applying the resulting control

%reporting task 2 again

w1 = 1;
w2 = 0.05;
w3 = 100;
N=10;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
x=[0.9492 0.43]';
Xsim2=zeros(2,Sim_time);
Xsim2(:,1)=x;
Usim2=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim2(:,i),N,xr1,w1,w2,w3,1);
        Usim2(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim2(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim2(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim2(i)=u;
        Xsim2(:,i+1)=x;
    end
end
time=toc
%applying noise: gaussian from a standard distribution multiplied by the
%standard  deviation required to have the same distribution than directly
%sampling from the wanted distribution

x=[0.9492 0.43]';
XsimN=zeros(2,Sim_time);
XsimN(:,1)=x;
UsimN=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        %creating the noise each step
        std_dev=sqrt(XsimN(2,i)^2/1000);
        noise=randn(1)*std_dev;

        u=NMPC(@CSTR,XsimN(:,i)+noise,N,xr1,w1,w2,w3,1);
        UsimN(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        XsimN(:,i+1)=x;
    else
        std_dev=sqrt(XsimN(2,i)^2/1000);
        noise=randn(1)*std_dev;

        u=NMPC(@CSTR,XsimN(:,i)+noise,N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        UsimN(i)=u;
        XsimN(:,i+1)=x;
    end
end
timeN=toc;
c=linspace(1,Sim_time,Sim_time)/5;
figure;
subplot(3, 1, 1);
plot(c,Xsim2(2,:),'r','DisplayName','Task2')
hold on;
plot(c,XsimN(2,:),'b','DisplayName','Task6')
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'m','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'m','LineStyle','--','DisplayName','x_r2')
legend;
xlabel('Time(s)');
ylabel('Concentration level');
grid on

subplot(3,1,2)
plot(c,Xsim2(1,:),'r','DisplayName','Task2')
hold on
plot(c,XsimN(1,:),'b','DisplayName','Task6')
legend;
xlabel('Time(s)');
ylabel('Temperature level');
grid on

subplot(3,1,3)
plot(c(1:end-1),Usim2,'r','DisplayName','Task2')
hold on
plot(c(1:end-1),UsimN,'b','DisplayName','Task6')
xlabel('Time (s)');
ylabel('u(k)');
legend;
grid on;


%redo everything from task 6 but reducing N
w1 = 1;
w2 = 0.05;
w3 = 100;
N=5;
Ts=0.2;
teta=20;
K=300;
M=5;
xf=0.3947;
xc=0.3816;
alpha=0.117;
xr1=0.6519;
xr2=0.5;
Sim_time=300/0.2;
x=[0.9492 0.43]';
Xsim2=zeros(2,Sim_time);
Xsim2(:,1)=x;
Usim2=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        u=NMPC(@CSTR,Xsim2(:,i),N,xr1,w1,w2,w3,1);
        Usim2(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Xsim2(:,i+1)=x;
    else
        u=NMPC(@CSTR,Xsim2(:,i),N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        Usim2(i)=u;
        Xsim2(:,i+1)=x;
    end
end
time=toc
%applying noise: gaussian from a standard distribution multiplied by the
%standard  deviation required to have the same distribution than directly
%sampling from the wanted distribution

x=[0.9492 0.43]';
XsimN=zeros(2,Sim_time);
XsimN(:,1)=x;
UsimN=zeros(1,Sim_time-1);
tic
for i=1:Sim_time-1
    disp(i)
    if i<=Sim_time/2
        %creating the noise each step
        std_dev=sqrt(XsimN(2,i)^2/1000);
        noise=randn(1)*std_dev;

        u=NMPC(@CSTR,XsimN(:,i)+noise,N,xr1,w1,w2,w3,1);
        UsimN(i)=u;
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        XsimN(:,i+1)=x;
    else
        std_dev=sqrt(XsimN(2,i)^2/1000);
        noise=randn(1)*std_dev;

        u=NMPC(@CSTR,XsimN(:,i)+noise,N,xr2,w1,w2,w3,1);
        x=CSTR(x,Ts,teta,K,xf,M,alpha,xc,u);
        UsimN(i)=u;
        XsimN(:,i+1)=x;
    end
end
timeN=toc;
c=linspace(1,Sim_time,Sim_time)/5;
figure;
subplot(3, 1, 1);
plot(c,Xsim2(2,:),'r','DisplayName','Task2')
hold on;
plot(c,XsimN(2,:),'b','DisplayName','Task6')
plot(c(1:Sim_time/2),ones(Sim_time/2,1)*xr1,'m','LineStyle','--','DisplayName','x_r1')
plot(c(Sim_time/2+1:end),ones(Sim_time/2,1)*xr2,'m','LineStyle','--','DisplayName','x_r2')
legend;
xlabel('Time(s)');
ylabel('Concentration level');
grid on

subplot(3,1,2)
plot(c,Xsim2(1,:),'r','DisplayName','Task2')
hold on
plot(c,XsimN(1,:),'b','DisplayName','Task6')
legend;
xlabel('Time(s)');
ylabel('Temperature level');
grid on

subplot(3,1,3)
plot(c(1:end-1),Usim2,'r','DisplayName','Task2')
hold on
plot(c(1:end-1),UsimN,'b','DisplayName','Task6')
xlabel('Time (s)');
ylabel('u(k)');
legend;
grid on;