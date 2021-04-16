%% Replicator Dynamics : ODE
clear, clc
b=1.9; Payoff_PD=[1 0; b 0]
dxdt =@(t,x) [x(1)^2*(1-x(1)-b*x(2)); x(1)*x(2)*(b*(1-x(2))-x(1))];
[t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
figure(1)
subplot(1,2,1)
plot(t,y(:,1),'r','LineWidth',2);grid
hold on
b=1.9; C=1; V=b; Payoff_HD=[(V-C)/2 V; 0 V/2]
dxdt =@(t,x) [x(1)*((V-C)/2*x(1)*(1-x(1))+V*x(2)*(1-x(1)-x(2)/2)); x(2)*(V/2*x(2)-(V-C)/2*x(1)^2-V*x(1)*x(2)-V/2*x(2)^2)];
[t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
plot(t,y(:,1),'m','LineWidth',2);
b=1.9; c=1; Payoff_SD=[(b-c)/2 b/2-c; b/2 0]
dxdt =@(t,x) [x(1)*((b-c)/2*x(1)*(1-x(1))+(c-b)*x(1)*x(2)+(b/2-c)*x(2)); x(2)*(b/2*x(1)*(1-x(1))+c/2*x(1)^2+(c-b)*x(1)*x(2))];
[t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
plot(t,y(:,1),'c','LineWidth',2);
title('Three Games in ODEs')
legend('Prisoners','Hawk Dove','Snowdrift')
xlabel('Time'); ylabel('Proportion of cooperators');
ylim([-0.05 1.05])
xlim([0 50])
hold off
%% Replicator Dynamics : PDE
clear, clc
rng('default')
% parameters
b=1.9; h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C=zeros(m,l,n);D=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C(:,:,1)=ones(m,l)*0.9+r*0.4;
D(:,:,1)=ones(m,l)*0.1+r*0.4;

% using finite difference scheme
for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C(1,k,i)=C(2,k,i); C(j,1,i)=C(j,2,i); C(1,1,i)=C(2,2,i); C(1,l,i)=C(2,l-1,i);
            C(m,k,i)=C(m-1,k,i); C(j,l,i)=C(j,l-1,i); C(m,l,i)=C(m-1,l-1,i); C(m,1,i)=C(m-1,2,i);
            D(1,k,i)=D(2,k,i); D(j,1,i)=D(j,2,i); D(1,1,i)=D(2,2,i); D(1,l,i)=D(2,l-1,i);
            D(m,k,i)=D(m-1,k,i); D(j,l,i)=D(j,l-1,i); D(m,l,i)=D(m-1,l-1,i); D(m,1,i)=D(m-1,2,i);
            
            % Cooperators
            C(j,k,i+1)=C(j,k,i)+dt*((C(j+1,k,i)+C(j-1,k,i)+C(j,k+1,i)+C(j,k-1,i)-4*C(j,k,i))/h^2)+...
                dt*C(j,k,i)^2*((1-C(j,k,i))-b*D(j,k,i));
            % Defectors
            D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                dt*C(j,k,i)*D(j,k,i)*(b*(1-D(j,k,i))-C(j,k,i));
        end
    end
end

% graphs
subplot(1,2,2)
vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
vD=sum(sum(D))/m^2; volD(:)=vD(1,1,:);
plot(t,volC,'r','LineWidth',2); grid
hold on

b=1.9; c=1; v=b; fr=(v-c)/2;
for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C(1,k,i)=C(2,k,i); C(j,1,i)=C(j,2,i); C(1,1,i)=C(2,2,i); C(1,l,i)=C(2,l-1,i);
            C(m,k,i)=C(m-1,k,i); C(j,l,i)=C(j,l-1,i); C(m,l,i)=C(m-1,l-1,i); C(m,1,i)=C(m-1,2,i);
            D(1,k,i)=D(2,k,i); D(j,1,i)=D(j,2,i); D(1,1,i)=D(2,2,i); D(1,l,i)=D(2,l-1,i);
            D(m,k,i)=D(m-1,k,i); D(j,l,i)=D(j,l-1,i); D(m,l,i)=D(m-1,l-1,i); D(m,1,i)=D(m-1,2,i);
            
            % Cooperators
            C(j,k,i+1)=C(j,k,i)+dt*((C(j+1,k,i)+C(j-1,k,i)+C(j,k+1,i)+C(j,k-1,i)-4*C(j,k,i))/h^2)+...
                dt*C(j,k,i)*(fr*C(j,k,i)*(1-C(j,k,i))+v*D(j,k,i)*(1-C(j,k,i)-D(j,k,i)/2));
            % Defectors
            D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                dt*D(j,k,i)*(v/2*D(j,k,i)-fr*C(j,k,i)^2-v*C(j,k,i)*D(j,k,i)-v/2*D(j,k,i)^2);
        end
    end
end

vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
vD=sum(sum(D))/m^2; volD(:)=vD(1,1,:);
plot(t,volC,'m','LineWidth',2);

b=1.9;  c=1;
for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C(1,k,i)=C(2,k,i); C(j,1,i)=C(j,2,i); C(1,1,i)=C(2,2,i); C(1,l,i)=C(2,l-1,i);
            C(m,k,i)=C(m-1,k,i); C(j,l,i)=C(j,l-1,i); C(m,l,i)=C(m-1,l-1,i); C(m,1,i)=C(m-1,2,i);
            D(1,k,i)=D(2,k,i); D(j,1,i)=D(j,2,i); D(1,1,i)=D(2,2,i); D(1,l,i)=D(2,l-1,i);
            D(m,k,i)=D(m-1,k,i); D(j,l,i)=D(j,l-1,i); D(m,l,i)=D(m-1,l-1,i); D(m,1,i)=D(m-1,2,i);
            
            % Cooperators
            C(j,k,i+1)=C(j,k,i)+dt*((C(j+1,k,i)+C(j-1,k,i)+C(j,k+1,i)+C(j,k-1,i)-4*C(j,k,i))/h^2)+...
                dt*C(j,k,i)*((b-c)/2*C(j,k,i)*(1-C(j,k,i))+(c-b)*C(j,k,i)*D(j,k,i)+(b/2-c)*D(j,k,i));
            % Defectors
            D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                dt*D(j,k,i)*(b/2*C(j,k,i)*(1-C(j,k,i))+c/2*C(j,k,i)^2+(c-b)*C(j,k,i)*D(j,k,i));
        end
    end
end

vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
vD=sum(sum(D))/m^2; volD(:)=vD(1,1,:);
plot(t,volC,'c','LineWidth',2);
title('Three Games in PDEs')
xlabel('Time'); ylabel('Proportion of cooperators');
legend('Prisoners','Hawk-Dove','Snowdrift')
ylim([-0.05 1.05])
xlim([0 50])
hold off

% figure(2)
% for i=1:n-1
%     gC(:,:)=C(:,:,i); gD(:,:)=D(:,:,i);
%     s1=surf(x,y,gC); s1.EdgeColor='g'; s1.FaceColor='k';
%     hold on
%     s2=surf(x,y,gD); s2.EdgeColor='b'; s2.FaceColor='k';
%     zlim([0 1]); title(sprintf('Time step = %d',i));
%     hold off
%     colormap bone
%     %     drawnow
%     pause
% end


%%
clear, clc
figure(2)
rng('default')
b=1.9; h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C=zeros(m,l,n);D=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C1(:,:,1)=ones(m,l)*0.9+r*0.4;
D1(:,:,1)=ones(m,l)*0.1+r*0.4;

for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C1(1,k,i)=C1(2,k,i); C1(j,1,i)=C1(j,2,i); C1(1,1,i)=C1(2,2,i); C1(1,l,i)=C1(2,l-1,i);
            C1(m,k,i)=C1(m-1,k,i); C1(j,l,i)=C1(j,l-1,i); C1(m,l,i)=C1(m-1,l-1,i); C1(m,1,i)=C1(m-1,2,i);
            D1(1,k,i)=D1(2,k,i); D1(j,1,i)=D1(j,2,i); D1(1,1,i)=D1(2,2,i); D1(1,l,i)=D1(2,l-1,i);
            D1(m,k,i)=D1(m-1,k,i); D1(j,l,i)=D1(j,l-1,i); D1(m,l,i)=D1(m-1,l-1,i); D1(m,1,i)=D1(m-1,2,i);
            
            % Cooperators
            C1(j,k,i+1)=C1(j,k,i)+dt*((C1(j+1,k,i)+C1(j-1,k,i)+C1(j,k+1,i)+C1(j,k-1,i)-4*C1(j,k,i))/h^2)+...
                dt*C1(j,k,i)^2*((1-C1(j,k,i))-b*D1(j,k,i));
            % Defectors
            D1(j,k,i+1)=D1(j,k,i)+dt*((D1(j+1,k,i)+D1(j-1,k,i)+D1(j,k+1,i)+D1(j,k-1,i)-4*D1(j,k,i))/h^2)+...
                dt*C1(j,k,i)*D1(j,k,i)*(b*(1-D1(j,k,i))-C1(j,k,i));
        end
    end
end
subplot(1,3,1)
gC1(:,:)=C1(:,:,46); gD1(:,:)=D1(:,:,46);
s1=surf(x,y,gC1); s1.EdgeColor='r'; s1.FaceColor='k';
hold on
s2=surf(x,y,gD1); s2.EdgeColor='b'; s2.FaceColor='k';
zlim([0 1]); title('Prisoners Dilemma at time = 3.6s');
xlabel('x-axis');ylabel('y-axis');zlabel('C and D');
hold off
colormap bone
% legend('Cooperator','Defector','Location','West')

%%

b=1.9; c=1; v=b; fr=(v-c)/2;
h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C2=zeros(m,l,n);D2=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C2(:,:,1)=ones(m,l)*0.9+r*0.4;
D2(:,:,1)=ones(m,l)*0.1+r*0.4;

for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C2(1,k,i)=C2(2,k,i); C2(j,1,i)=C2(j,2,i); C2(1,1,i)=C2(2,2,i); C2(1,l,i)=C2(2,l-1,i);
            C2(m,k,i)=C2(m-1,k,i); C2(j,l,i)=C2(j,l-1,i); C2(m,l,i)=C2(m-1,l-1,i); C2(m,1,i)=C2(m-1,2,i);
            D2(1,k,i)=D2(2,k,i); D2(j,1,i)=D2(j,2,i); D2(1,1,i)=D2(2,2,i); D2(1,l,i)=D2(2,l-1,i);
            D2(m,k,i)=D2(m-1,k,i); D2(j,l,i)=D2(j,l-1,i); D2(m,l,i)=D2(m-1,l-1,i); D2(m,1,i)=D2(m-1,2,i);
            
            % Cooperators
            C2(j,k,i+1)=C2(j,k,i)+dt*((C2(j+1,k,i)+C2(j-1,k,i)+C2(j,k+1,i)+C2(j,k-1,i)-4*C2(j,k,i))/h^2)+...
                dt*C2(j,k,i)*(fr*C2(j,k,i)*(1-C2(j,k,i))+v*D2(j,k,i)*(1-C2(j,k,i)-D2(j,k,i)/2));
            % Defectors
            D2(j,k,i+1)=D2(j,k,i)+dt*((D2(j+1,k,i)+D2(j-1,k,i)+D2(j,k+1,i)+D2(j,k-1,i)-4*D2(j,k,i))/h^2)+...
                dt*D2(j,k,i)*(v/2*D2(j,k,i)-fr*C2(j,k,i)^2-v*C2(j,k,i)*D2(j,k,i)-v/2*D2(j,k,i)^2);
        end
    end
end

subplot(1,3,2)
    gC2(:,:)=C2(:,:,46); gD2(:,:)=D2(:,:,46);
    s1=surf(x,y,gC2); s1.EdgeColor='m'; s1.FaceColor='k';
    hold on
    s2=surf(x,y,gD2); s2.EdgeColor='b'; s2.FaceColor='k';
    zlim([0 1]); title('Hawk-Dove at time = 3.6s');
    xlabel('x-axis');ylabel('y-axis');zlabel('C and D');
    hold off
    colormap bone

%%
b=1.9;  c=1;
h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C3=zeros(m,l,n);D3=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C3(:,:,1)=ones(m,l)*0.9+r*0.4;
D3(:,:,1)=ones(m,l)*0.1+r*0.4;

for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % periodic boundary conditions
            C3(1,k,i)=C3(2,k,i); C3(j,1,i)=C3(j,2,i); C3(1,1,i)=C3(2,2,i); C3(1,l,i)=C3(2,l-1,i);
            C3(m,k,i)=C3(m-1,k,i); C3(j,l,i)=C3(j,l-1,i); C3(m,l,i)=C3(m-1,l-1,i); C3(m,1,i)=C3(m-1,2,i);
            D3(1,k,i)=D3(2,k,i); D3(j,1,i)=D3(j,2,i); D3(1,1,i)=D3(2,2,i); D3(1,l,i)=D3(2,l-1,i);
            D3(m,k,i)=D3(m-1,k,i); D3(j,l,i)=D3(j,l-1,i); D3(m,l,i)=D3(m-1,l-1,i); D3(m,1,i)=D3(m-1,2,i);
            
            % Cooperators
            C3(j,k,i+1)=C3(j,k,i)+dt*((C3(j+1,k,i)+C3(j-1,k,i)+C3(j,k+1,i)+C3(j,k-1,i)-4*C3(j,k,i))/h^2)+...
                dt*C3(j,k,i)*((b-c)/2*C3(j,k,i)*(1-C3(j,k,i))+(c-b)*C3(j,k,i)*D3(j,k,i)+(b/2-c)*D3(j,k,i));
            % Defectors
            D3(j,k,i+1)=D3(j,k,i)+dt*((D3(j+1,k,i)+D3(j-1,k,i)+D3(j,k+1,i)+D3(j,k-1,i)-4*D3(j,k,i))/h^2)+...
                dt*D3(j,k,i)*(b/2*C3(j,k,i)*(1-C3(j,k,i))+c/2*C3(j,k,i)^2+(c-b)*C3(j,k,i)*D3(j,k,i));
        end
    end
end

subplot(1,3,3)
    gC3(:,:)=C3(:,:,46); gD3(:,:)=D3(:,:,46);
    s1=surf(x,y,gC3); s1.EdgeColor='c'; s1.FaceColor='k';
    hold on
    s2=surf(x,y,gD3); s2.EdgeColor='b'; s2.FaceColor='k';
    zlim([0 1]); title('Snowdrift at time = 3.6s');
    xlabel('x-axis');ylabel('y-axis');zlabel('C and D');
    hold off
    colormap bone

