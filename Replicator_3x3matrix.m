clc, clear

% Initial setting
x0=[0.1 0.3 0.6]; h=0.005; tspan=[0 40];

% Numerical approximation for RPS dynamics
[t,y]=ode45(@RPS,tspan,x0,h);

% Plot for replicator dynamics with unit simplex
figure(2)
subplot(1,2,1)
plot3(y(:,1),y(:,2),y(:,3),'b')
title('RPS-game in Phase Plane with Simplex')
hold on; 
vertex=[1 0 0; 0 1 0; 0 0 1]; face=[1 2 3];
patch('Faces',face,...
    'Vertices',vertex,'Facecolor',[0.98 0.98 0.98]);
view(-45,0)
xticklabels({'','','(Scissors)'})
yticklabels({'','','(Rock)'})
zticklabels({'','','','','','','','','','','(Paper)'})
hold off; axis off

% Plot for the proportion of the cooperation
subplot(1,2,2)
plot(t,y(:,1),'r',t,y(:,2),'b',t,y(:,3),'k')
legend('Rock','Paper','Scissors')
title('Comparison of Three Different Strategies')
xlabel('Time (iteration)')
ylabel('Proportion of the cooperation')
grid on
