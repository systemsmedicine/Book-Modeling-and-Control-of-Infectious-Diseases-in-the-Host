%% Modeling and Control of Infectious Diseases in the Host with MATLAB and R
%% 6.10 SIMULATING HIV MODEL
%% Page 123

clc; clear all; close all;
Tmax=10*365;

%% Solving Differential Equations
[t,y]=ode45(@models,[0 Tmax],[1000 0 150 0 1e-3]);% ODE function
T=y(:,1);    % Healthy (uninfected) cells 
Ti=y(:,2);   % Infected cells (Active)
M=y(:,3);    % Macrophages
Mi=y(:,4);   % Infected Macrophages
V=y(:,5);    % Virus

% Plotting
figure(1); clf;
sim_V=plot(t/365,T,'LineWidth',2);
axis([0,Tmax/365,1,1000]);
xlabel('Time (days)');
ylabel('T(cells/mm^3)','fontsize',10.0)
hold on
figure(2); clf;
plot(t/365,V/0.001,'LineWidth',2);ylabel('V (copies/mm^3)')
xlabel('Time (years)');
axis([0,Tmax/365,1,2500000]);
figure(3); clf;
plot(t/365,Ti,'LineWidth',2);ylabel('T^* (cells/mm^3)')
xlabel('Time (years)');
figure(4); clf;
plot(t/365,M,t/365,Mi,'LineWidth',2); ylabel('M (cells/mm^3)')
xlabel('Time (years)');

