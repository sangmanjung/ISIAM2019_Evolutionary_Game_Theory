%% Replicator Dynamics : ODE
clear, clc
% c is a costs, b is a benefits
b=1.9; c=1; Payoff=[(b-c)/2 b/2-c; b/2 0]
dxdt =@(t,x) [x(1)*((b-c)/2*x(1)*(1-x(1))+(c-b)*x(1)*x(2)+(b/2-c)*x(2)); x(2)*(b/2*x(1)*(1-x(1))+c/2*x(1)^2+(c-b)*x(1)*x(2))];
[t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
figure(1)
subplot(1,2,1)
plot(t,y(:,1),'g',t,y(:,2),'b','LineWidth',2);grid
legend('Cooperators','Defectors')
title('Snowdrift games in ODE')
xlabel('Time'); ylabel('C and D');
%% Replicator Dynamics : PDE
clear, clc
figure(1)
rng('default')
% parameters
b=1.9;  c=1; h=2/sqrt(10); dt=(h^2)/5;
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
                dt*C(j,k,i)*((b-c)/2*C(j,k,i)*(1-C(j,k,i))+(c-b)*C(j,k,i)*D(j,k,i)+(b/2-c)*D(j,k,i));
            % Defectors
            D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                dt*D(j,k,i)*(b/2*C(j,k,i)*(1-C(j,k,i))+c/2*C(j,k,i)^2+(c-b)*C(j,k,i)*D(j,k,i));
        end
    end
end

% graphs
subplot(1,2,2)
vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
vD=sum(sum(D))/m^2; volD(:)=vD(1,1,:);
plot(t,volC,'g',t,volD,'b','LineWidth',2); grid
legend('Cooperator','Defector')
title('Snowdrift games in PDE'); xlabel('Time'); ylabel('C and D');
% ylim([0 1])

% figure(2)
% for i=1:n-1
%     gC(:,:)=C(:,:,i); gD(:,:)=D(:,:,i);
%     s1=surf(x,y,gC); s1.EdgeColor='g'; s1.FaceColor='k';
%     hold on
%     s2=surf(x,y,gD); s2.EdgeColor='b'; s2.FaceColor='k';
%     zlim([0 1]); title(sprintf('Time step = %d',i));
%     hold off
%     colormap bone
%         drawnow
% %     pause
% end
