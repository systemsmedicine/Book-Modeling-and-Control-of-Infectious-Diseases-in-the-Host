%% Modeling and Control of Infectious Diseases in the Host with MATLAB and R
%% 8.13 MATLAB CODE FOR ALGORITHM 3
%% Page 188

clear all; close all; clc;
tic  % Start Simulation Time
n=4; % State dimension
m=2; % Number of matrices
cost=[1,1,1,1]'; % Cost vector
T_max=336; % Final Time
dTsim=28;  % Simulation Time
%% Steps
N=T_max/dTsim; t=linspace(0,T_max,N+1);
%% Set up Model parameters.
mu=1e-04;                            % Mutation rate
M=[0 1 1 0;1 0 0 1;1 0 0 1;0 1 1 0]; % Mutation pattern
death_rate=0.24;
%% Reproduction rates
rho=[0.05, 0.05; 0.27, 0.15;0.01,0.25;0.27,0.27];
R1=diag(rho(:,1));  R2=diag(rho(:,2));
%% State matrices
beta=0;
A(:,:,1)=expm((R1-death_rate*eye(size(R1))+mu*M-beta*eye(size(R1)))*dTsim);
A(:,:,2)=expm((R2-death_rate*eye(size(R2))+mu*M-beta*eye(size(R1)))*dTsim);
%% Initial Conditions for the State
Vinit=1000; xt1=Vinit; xt2=xt1*mu; xt3=xt1*mu; xt4=xt2*mu+xt3*mu;
x0=[xt1 xt2 xt3 xt4]';
%% Algorithm Initialization
Psi=x0; sequence=0;
lb = zeros(length(A),1); % Lower Bound of A
T_LP=1e-5;   % Tolerance for Linear Programming
T_kc=1+1e-7; % Tolerance to keep columns
% Keep count of how many linprog warnings
maxiterr_count=0;
linprogerr_count=0;
%% COMPUTING CLEANED COLUMNS FORWARD FOR N/2 STEPS
for i=1:(N/2)
   % Generating Psi
   [dummy,dim]=size(Psi); LastPsi=Psi(:,1:dim);
   Psi=[A(:,:,1)'*LastPsi A(:,:,2)'*LastPsi]; 
   % Following The Sequences 
   clear c1 c2  
   for r=1:dim
        c1(r,:)=[1,sequence(r,:)];
        c2(r,:)=[2,sequence(r,:)];
   end
   sequence=[c1;c2];
   % Linear Programming Inequality
    Ain=-Psi';
   [raws_A,columns_a]=size(Ain);
    % Computing myu
    p=0;
    clear myu;
    delate_columns=0; jf=0;
    %% Cycle for clear columns
    for j=1:raws_A
        % Following which column is for the new Matrix
        jp=j-delate_columns;
        % Delating column just for Optimization
        Ak=Ain;  Ak(jp,:)=[];
        [raws_Ak,columns_Ak]=size(Ak);
        b=-ones(raws_Ak,1);
        % Linear Programming Function
        options = optimset('TolFun',T_LP);
        [xmin(:,j),fval,exitflag,output,lambda] = linprog(Psi(:,j)',Ak,b,[],[],lb,[],[],options);
        if exitflag == 0
           maxiterr_count=maxiterr_count+1;
        end
        if exitflag < 0
            linprogerr_count=linprogerr_count+1;
        end
        % Miu Computation
        myu(j)=Psi(:,j)'*xmin(:,j);
        % Cleaning Columns for the next step
        if ((myu(j))<T_kc)
           p=p+1;
           Psiclean(:,p)=Psi(:,j);
           cleansequence(p,:)=sequence(j,:);
        end
        % Delating Columns for LP
        if ((myu(j))>=T_kc)
           Ain(jp,:)=[];
           delate_columns=delate_columns+1;
        end
     end
    % Psi and Sequence for Next Step
    Psi_before=Psi;
    sequence_before=sequence;
    clear Psi sequence 
    Psi=Psiclean;
    sequence=cleansequence;
    clear cleansequence Psiclean LastPsi 
    %% Print how many columns are at the end of the optimization
    [rf,cf]=size(Psi);
    fprintf(' The simulation is in the step %i and the columns are %i\n',i,cf)
end
%% COMPUTING CLEANED COLUMNS BACKWARD FOR N/2 STEPS  
Psik=cost; sequencek=0;
for ik=1:(N/2)
   % Generating Psi
   [dummyk,dimk]=size(Psik);
   LastPsik=Psik(:,1:dimk);
   Psik=[A(:,:,1)'*LastPsik A(:,:,2)'*LastPsik]; 
   % Following The Sequences 
   clear ck1 ck2  
   for rk=1:dimk
        ck1(rk,:)=[1,sequencek(rk,:)];
        ck2(rk,:)=[2,sequencek(rk,:)];
   end
   sequencek=[ck1;ck2];
   % Linear Programming Inequality
   Aink=-Psik';
   [raws_Ak,columns_ak]=size(Aink);
   % Computing myu
    pk=0;
    clear myuk;
    delate_columnsk=0;
    jfk=0;
   % Cycle for clear columns
    for jk=1:raws_Ak    
        % Following which column is for the new Matrix
        jpk=jk-delate_columnsk;
        % Delating column just for Optimization
        Akk=Aink;
        Akk(jpk,:)=[];
        [raws_Akk,columns_Akk]=size(Akk);
        bk=-ones(raws_Akk,1);
        % Linear Programming Function
        options = optimset('TolFun',T_LP);
        [xmink(:,jk),fval,exitflag,output,lambda] = linprog(Psik(:,jk)',Akk,bk,[],[],lb,[],[],options);
        % Miu Computation
        myuk(jk)=Psik(:,jk)'*xmink(:,jk);
        % Cleaning Columns for the next step
        if ((myuk(jk))<T_kc)
           pk=pk+1;
           Psicleank(:,pk)=Psik(:,jk);
           cleansequencek(pk,:)=sequencek(jk,:);
        end
        % Delating Columns for LP
        if ((myuk(jk))>=T_kc)
           Aink(jpk,:)=[];
           delate_columnsk=delate_columnsk+1;
        end
    end
    %% Psi and Sequence for Next Step
    Psi_beforek=Psik;
    sequence_beforek=sequencek;
    clear Psik sequencek
    Psik=Psicleank;
    sequencek=cleansequencek;
    clear cleansequencek Psicleank LastPsik 
    %% Print how many columns are at the end of the optimization
    [rfk,cfk]=size(Psik);
    fprintf(' The simulation is in the step %i and the columns are %i\n',ik,cfk)

end
%% COUPLED FORWARD AND BACKWARDS SEQUENCES 
[r_opt,c_opt]= size(Psi); [r_optk,c_optk]= size(Psik);
gp=0;
for i=1:c_opt
    for j=1:c_optk
        gp=gp+1;
        cost_final(gp)=Psi(:,i)'*Psik(:,j);
        eye_sequence(gp,:)=[i,j];
    end
end
%% Looking the Optimal Coupled Sequence
[val_opt,seq_opt]=min(cost_final);
%% Optimal Sequence
opt_seq=eye_sequence(seq_opt,:);
optimal_control=[fliplr(sequence(opt_seq(1),1:(N/2))), sequencek(opt_seq(2),1:(N/2))];
clc;
%% Printing Simulation Time in Seconds
display('-------------------------------------------------------------- ');
display('Simulation Time'); toc
display('-------------------------------------------------------------- ');
%% Printing the total number of columns
fprintf('  The total of clean columns in the Forward is %i \n ',c_opt);
fprintf(' The total of clean columns in the Backward is %i\n ',c_optk);
fprintf(' The total of clean columns for %i stetps is %i\n ',N,(c_opt*c_optk));
display('--------------------------------------------------------------');
%% Printing the total Viral Load in the last Step
fprintf('  The total Viral Load in the step %i is %f \n ',N, val_opt);
display('--------------------------------------------------------------');
%% Print some things to do with errors
fprintf('No. of lin prog errors is %i \n',linprogerr_count);
fprintf('No. of lin prog max iteration warnings is %i \n',maxiterr_count);
%% Plotting Optimal Control
figure(1); clf;
stairs(t(:,1:N),optimal_control); 
ylabel('\sigma'); xlabel('Time (days)'); axis([0,T_max,0,3]);
title('Optimal Control');

