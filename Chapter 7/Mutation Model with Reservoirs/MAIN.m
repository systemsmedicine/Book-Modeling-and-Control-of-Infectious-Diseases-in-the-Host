%% Modeling and Control of Infectious Diseases in the Host with MATLAB and R
%% 7.12 SIMULATING THE MUTATION MODEL (7.5)
%% Page 144

clear all; close all; clc

%% State matrices
[A1,C1]=build_linear_model(1);
[A2,C2]=build_linear_model(2);

%% Initial Viral Load
nx=length(C1);
Vinit=600;
x0=zeros(nx,1);
x0(1)=150;
x0(2)=50;
x0(3)=Vinit;

%% Simulation time
dTsim=30;              % Time step for simulations
dTc=30;                % Time step for decisions
Tmax=12*30*16;         % Max simulation time
t=0:dTsim:Tmax; t=t';  % Time vector
q=C1';                 % Cost vector
N=Tmax/dTc;
%% Initialise data variables
nt=length(t); xt=zeros(nt,nx); w=zeros(nt,1); c=w; xt(1,:)=x0'; z=xt(1,:)';
Phi1=expm(A1*dTsim);
Phi2=expm(A2*dTsim);
failed=0;
c(1)=1;
g=1;
t_lastc=0; % time at which the last control was computed
if g==1 Phi=Phi1; else Phi=Phi2; end

%% Simulate Viral Escape
for k=1:length(t)-1;
    z=Phi*z;
    xt(k+1,:)=z';
    if (t(k+1)-t_lastc) >= dTc % time to compute a new decision
        t_lastc=t(k+1);
        [g,failed]=vir_escape(q'*z,failed,1000);
        if g==1
            Phi=Phi1;
        else
            Phi=Phi2;
        end
    end
    c(k+1)=g;
end

%% Record data for display
C_Ti=zeros(1,27); for i=0:8; C_Ti(3*i+1)=1; end % Infected Cells
C_L=zeros(1,27); for i=0:8; C_L(3*i+2)=1; end   % Latently Infected Cells
v_vir=xt*q; c_vir=c;
Ti=xt*C_Ti';
TL=xt*C_L';
c_vir=c;    % Switched Rule
Result_escape= v_vir(length(v_vir));

%% Figure
figure(1); clf;
subplot(3,1,1);
semilogy(t/(30*12),v_vir);
ylabel('copies/ml');
subplot(3,1,2);
semilogy(t/(30*12),Ti,t/(30*12),TL );
legend('T_i','T_L')
ylabel('cells/mm^3');
subplot(3,1,3);
stairs(t/(30*12),c_vir); ylabel('\sigma(k)');
axis([0,Tmax/(30*12),0,3])
xlabel('Time (years)');


