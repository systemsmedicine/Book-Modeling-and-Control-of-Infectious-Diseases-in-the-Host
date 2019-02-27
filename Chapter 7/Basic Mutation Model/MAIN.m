%% Modeling and Control of Infectious Diseases in the Host with MATLAB and R
%% 7.11 SIMULATING THE MUTATION MODEL (7.2) FOR SCENARIO 1
%% Page 141

clear all; close all; clc
%% Set up system parameters.
mu=1e-04;                            % Mutation rate
M=[0 1 1 0;1 0 0 1;1 0 0 1;0 1 1 0]; % Mutation pattern
death_rate=0.24;                     % Viral clearance
Vinit=1000;                          % Initial condition

%% Time Simulation
dTsim=28;                            % Time step for simulations
Tmax=336;                            % Max simulation time
t=0:dTsim:Tmax; t=t';                % Time vector

%% Reproduction rates: Scenario 1
rho=[0.05,0.05;0.27,0.05;0.05, 0.27;0.27,0.27];    
R1=diag(rho(:,1)); R2=diag(rho(:,2));

%% State matrices
A1=R1-death_rate*eye(size(R1))+mu*M;
A2=R2-death_rate*eye(size(R2))+mu*M;

%% Initial Conditions
x11(1)=Vinit; x21(1)=x11(1)*mu; x31(1)=x11(1)*mu; 
x41(1)=x21(1)*mu+x31(1)*mu;
z=[x11(1);x21(1);x31(1);x41(1)];
[s,g]=min([x21(1) x31(1)]);
c1(1)=g;
if g==1 At=A1; else At=A2; end

%% Main cycle for Virologic Failure
%% If you want to use Virological Failure Treatment uncomment lines 39-42 and comment line 44
for k=1:length(t)-1;
    z=expm(At*dTsim)*z;
    x11(k+1)=z(1);x21(k+1)=z(2);x31(k+1)=z(3);x41(k+1)=z(4);
    %% Virological Failure   
%     if sum(z) > 1000 
%         if c1(1)==1 g=2; end
%         if c1(1)==2 g=1; end
%     end 
    %% Virological Failure 
    if g==1 g=2; else g=1; end 
    c1(k+1)=g;
    if g==1 At=A1; else At=A2;end
end

%% Plotting
figure(1); clf;
subplot(2,1,1);
semilogy(t/28,x11,t/28,x21,t/28,x31,t/28,x41,t/28,x11+x21+x31+x41); 
legend('G1','G2','G3','G4','viral load');
xlabel('Time (months)'); ylabel('copies/ml');subplot(2,1,2);
stairs(t/28,c1);
ylabel('\sigma'); xlabel('Time (months)');axis([0,Tmax/28,0,3]);

