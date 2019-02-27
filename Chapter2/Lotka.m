%% Modeling and Control of Infectious Diseases in the Host with MATLAB and R
%% 2.7 SIMULATING LOTKA–VOLTERRA MODEL
%% Page 30

% Clearing window and memory
clc; clear all; close all;
% Simulation time
t0=200; 
% Initial conditions
x_0=10;   % Pray
z_0=1;    % Predator
%Solving odes
[t,y]=ode45(@funodes,[0 t0],[x_0 z_0]);
% changing variables for presentation
x=y(:,1); z=y(:,2);
% Plotting Results 
figure(1); 
plot(t,x,t,z,'LineWidth',2);
xlabel('Time','fontsize',20);
ylabel('Population number','fontsize',20);
legend({'Pray','Predator'},'FontSize',20);
ha1=gca; set(ha1,'LineWidth',2,'FontSize',20);
hold on
% The differential equations needs to write as a function either at
% the end of the script or in another file.

% REMARK: In MATLAB version R2016b and after you can have function definitions in a script

function  dy = funodes(t,y)
% Changing of variables for presentation
x=y(1);  
z=y(2);
% Model parameters
k1=1; k2=0.2; k3=0.05; k4=0.1;
% ODEs
dy = zeros(length(y),1);
dy(1) = k1*x -k2*x*z;
dy(2) = k3*x*z-k4*z;
end



