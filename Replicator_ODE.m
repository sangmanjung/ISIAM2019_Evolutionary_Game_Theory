%% Replicator Dynamics : ODE 
clear, clc
% parameters
h=0.1; t=0:h:50; n=length(t); b=1.9; C=b; V=1; a=1;
% information of the analytic expressions of replicator dynamics
disp('Prisoners Dilemma games :  dxdt=-x^2(1-x)(1-b), b is parameter')
disp('Bach or Stravinsky games : dxdt=-x(5x-2)(x-1)')
disp('Stag-Hunt games :          dxdt=-x(2x-1)(x-1)')
disp('Chicken games :            dxdt=x(5x-4)(x-1)')
disp('Hawk-Dove games :          dxdt=x/2(Cx(x-1)-V(x+1)), C and V are parameters')
disp('generalized CK  :          dxdt=x(1-x)(b(1-x)-a), a and b are parameters')
% New Hawk-Dove : dxdt=x(1-x)(Cx/2+v/2)
  
% information of the payoff matrix
PD=[1 0; b 0]
BS=[3 0; 0 2]
SH=[2 0; 1 1]
CK=[0 -1; 1 -5]
HD=[(V-C)/2 V; 0 V/2]
gen_CK=[0 -a; a -b]

% Graphs
figure(1)
% the frequency of strategy
y1(1)=0.3; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,1)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.3')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off

% the frequency of strategy
y1(1)=0.4; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,2)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.4')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off

% the frequency of strategy
y1(1)=0.41; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,3)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.41')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off

% the frequency of strategy
y1(1)=0.5; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,4)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.5')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off

% the frequency of strategy
y1(1)=0.51; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,5)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.51')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off

% the frequency of strategy
y1(1)=0.9; y2(1)=y1(1); y3(1)=y2(1); y4(1)=y3(1); y5(1)=y4(1); y6(1)=y5(1);
% numerical approximation (Euler's method)
for i=1:n-1
    y1(i+1)=y1(i)+h*(y1(i)^2*(1-y1(i))*(1-b)); % Prisoner's dilemma.
    y2(i+1)=y2(i)-h*(y2(i)*(5*y2(i)-2)*(y2(i)-1)); % Bach or Stravinsky
    y3(i+1)=y3(i)-h*(y3(i)*(y3(i)-1)*(2*y3(i)-1)); % Stag-Hunt
    y4(i+1)=y4(i)-h*(y4(i)*(y4(i)-1)*(5*y4(i)-4)); % Chicken
    y5(i+1)=y5(i)+h*(y5(i)/2*(C*y5(i)*(y5(i)-1)-V*(y5(i)+1))); % Hawk-Dove
    y6(i+1)=y6(i)+h*(y6(i)*(1-y6(i))*(b*(1-y6(i))-a)); % generalized CK
end
subplot(2,3,6)
plot(t,y1,'r','LineWidth',2.2)
hold on
plot(t,y2,'b','LineWidth',2.2)
plot(t,y3,'g','LineWidth',2.2)
plot(t,y4,'k','LineWidth',2.2)
plot(t,y5,'m','LineWidth',2.2)
plot(t,y6,'LineWidth',2.2)
legend('Prisoners','B or S','Stag-hunt','Chicken','Hawk-Dove','generalized CK')
title('Initial value : 0.9')
xlabel('Time (step size is 0.1)')
ylabel('Proportion of the cooperation')
ylim([-0.2 1.1])
grid on
hold off
 