function [A,C]=build_linear_model(sigma)
% 'sigma' = 1 or 2 denotes the drug combination used.
mu=3e-05;
%% Load some model parameters used in the calculations.
fitness1=[1;0.95;0.95]; fitness2=[1;0.95;0.95]; % Fitness 
% Drug efficiencies for each regimen
infection_eff1=[0.2;0.9;1]; production_eff1=[0.25;0.5;1];
infection_eff2=[0.2;0.5;1]; production_eff2=[0.1;0.8;1];
%% Initialazing
A=zeros(27,27); C=zeros(1,27);
for i=0:8; C(3*i+3)=1; end
def_params.KT=3.4714e-5;  def_params.dT=0.4;      def_params.PT=44;
def_params.KM=4.533e-7;   def_params.dM=0.001;    def_params.PM=44;
def_params.T=600;         def_params.M=140;       def_params.dV=2.4;
%% Build individual SS matrices
for i1=1:3
    for i2=1:3
        i=i1+3*i2-3; %Model subindex, i=1,2..9
        params=def_params;
        if sigma ==1
            params.KT = params.KT * infection_eff1(i1);
            params.KM = params.KM * infection_eff1(i1);
            params.PT = params.PT * fitness1(i1)*fitness2(i2)*production_eff1(i1);
            params.PM = params.PM * fitness1(i1)*fitness2(i2)*production_eff1(i1);
        else
            params.KT = params.KT * infection_eff2(i2);
            params.KM = params.KM * infection_eff2(i2);
            params.PT = params.PT * fitness1(i1)*fitness2(i2)*production_eff2(i2);
            params.PM = params.PM * fitness1(i1)*fitness2(i2)*production_eff2(i2);
        end
        [Ai,Bi,Ci]=mutant_model(params);
        CC(:,:,i)=Ci;
        A(i*3-2:i*3,i*3-2:i*3)=Ai;
    end
end
%% Define mutation matrix. First define edges in the graph
edges=[     1,2;    1,4;  2,3;    2,5;    3,6;
            4,5;    4,7;  5,6;    5,8;    6,9;
            7,8;    8,9];
for i=1:length(edges)
    i1=edges(i,1); i2=edges(i,2);
    N1=i1*3-2:i1*3; N2=i2*3-2:i2*3;
    A(N1,N2)=mu*CC(:,:,i2); A(N2,N1)=mu*CC(:,:,i1);
end
end
function [A,B,C]=mutant_model(params)
B=[1;0;0]; 
A=[
    -params.dT,      0,            params.KT*params.T;
    0,           -params.dM,       params.KM*params.M;
    params.PT,    params.PM,       -params.dV,];
C=A; C(1:3,1:2)=0; C(3,3)=0; 
end