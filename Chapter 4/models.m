function  dy= models(t,y)
s1=10;s2=0.15;
K2=4.5714e-005;K4=4.3333e-008;
K5=38;K6=35;
d1=0.01;d2=0.4;d3=0.001;d4=0.001;d5=2.4;
pT=0.01*(y(1)*y(5))/(y(5)+300);
pM=0.003*(y(3)*y(5))/(y(5)+220);
dy = zeros(5,1);
dy(1) = s1 -d1*y(1)-K2*y(5)*y(1)+pT;
dy(2) = K2*y(5)*y(1)-d2*y(2);
dy(3) = s2 -K4*y(3)*y(5)-d3*y(3)+pM;
dy(4) = K4*y(5)*y(3)-d4*y(4);
dy(5) = K5*y(2)+ K6*y(4)-d5*y(5);
end