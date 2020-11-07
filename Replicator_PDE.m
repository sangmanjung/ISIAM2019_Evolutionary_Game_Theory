%% Replicator Dynamics : PDE 
clear, clc
rng('default') 
% parameters
a=0; b=1.9;  C=b; V=1; h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:50; x=0:h:10; y=0:h:10;
n=length(t); m=length(x); l=length(y);
u1=zeros(m,l,n);u2=zeros(m,l,n);u3=zeros(m,l,n);
u4=zeros(m,l,n);u5=zeros(m,l,n);u6=zeros(m,l,n);

% payoffs
PD=[1 0; b 0]
BS=[3 0; 0 2]
SH=[2 0; 1 1]
CK=[0 -1; 1 -5]
HD=[(V-C)/2 V; 0 V/2]
gen_CK=[0 -a; a -b]

% set the initial condition
% for j=1:m
%     for k=1:l
%         if (abs(x(j))<=1) && (abs(y(k)) <=1) % later abs(x*y) cases
%                         u0(j,k)=cos(pi*x(j)*y(k)/2);
%             u0(j,k)=cos(pi*x(j)/2)*cos(pi*y(k)/2);
%         else
%             u0(j,k)=0;
%         end
%     end
% end
% 
% u1(:,:,1)=u0; u2(:,:,1)=u0; u3(:,:,1)=u0;
% u4(:,:,1)=u0; u5(:,:,1)=u0; u6(:,:,1)=u0;

u1(:,:,1)=rand(m,l); u2(:,:,1)=rand(m,l); u3(:,:,1)=rand(m,l);
u4(:,:,1)=rand(m,l); u5(:,:,1)=rand(m,l); u6(:,:,1)=rand(m,l);

% u1(:,:,1)=ones(m,l)*0.5; u2(:,:,1)=ones(m,l)*0.5;
% u3(:,:,1)=ones(m,l)*0.5; u4(:,:,1)=ones(m,l)*0.5;
% u5(:,:,1)=ones(m,l)*0.5; u6(:,:,1)=ones(m,l)*0.5;

% using finite difference scheme
for i=1:n-1
    for j=2:m-1
        for k=2:l-1
            % zero-flux boundary condition
            %             u(i,j,1)=u(i,j,3); u(i,j,l)=u(i,j,l-2);
            %             u(i,m,k)=u(i,m-2,k); u(i,1,k)=u(i,2,k);
            
            % FTCS scheme
            u1(j,k,i+1)=u1(j,k,i)+dt*((u1(j+1,k,i)+u1(j-1,k,i)+u1(j,k+1,i)+u1(j,k-1,i)-4*u1(j,k,i))/h^2)-...
                dt*u1(j,k,i)^2*(1-u1(j,k,i))*(1-b); % PD
            u2(j,k,i+1)=u2(j,k,i)+dt*((u2(j+1,k,i)+u2(j-1,k,i)+u2(j,k+1,i)+u2(j,k-1,i)-4*u2(j,k,i))/h^2)-...
                dt*u2(j,k,i)*(5*u2(j,k,i)-2)*(u2(j,k,i)-1); % B or S
            u3(j,k,i+1)=u3(j,k,i)+dt*((u3(j+1,k,i)+u3(j-1,k,i)+u3(j,k+1,i)+u3(j,k-1,i)-4*u3(j,k,i))/h^2)+...
                dt*u3(j,k,i)*(2*u3(j,k,i)-1)*(u3(j,k,i)-1); % Chicken
            u4(j,k,i+1)=u4(j,k,i)+dt*((u4(j+1,k,i)+u4(j-1,k,i)+u4(j,k+1,i)+u4(j,k-1,i)-4*u4(j,k,i))/h^2)-...
                dt*u4(j,k,i)*(2*u4(j,k,i)-1)*(u4(j,k,i)-1); % SH
            u5(j,k,i+1)=u5(j,k,i)+dt*((u5(j+1,k,i)+u5(j-1,k,i)+u5(j,k+1,i)+u5(j,k-1,i)-4*u5(j,k,i))/h^2)+...
                dt*u5(j,k,i)*(1-u5(j,k,i))*(C*u5(j,k,i)/2+V/2); % HD
            u6(j,k,i+1)=u6(j,k,i)+dt*((u6(j+1,k,i)+u6(j-1,k,i)+u6(j,k+1,i)+u6(j,k-1,i)-4*u6(j,k,i))/h^2)+...
                dt*u6(j,k,i)*(1-u6(j,k,i))*(b*(1-u6(j,k,i))-a); % generalized CK
        end
    end
end

% mesh graph
% figure(1)
% set(gcf, 'Position',  [350, 350, 500, 310])
% for i=1:n
%     gu1(:,:)=u1(:,:,i);
%     gu2(:,:)=u2(:,:,i);
%     gu3(:,:)=u3(:,:,i);
%     gu4(:,:)=u4(:,:,i);
%     gu5(:,:)=u5(:,:,i);
%     subplot(1,5,1)
%     mesh(x,y,gu1)
%     axis([-4 4 -4 4 0 1])
%     title(sprintf('Time step = %d',i));
%     xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
%     subplot(1,5,2)
%     mesh(x,y,gu2)
%     axis([-4 4 -4 4 0 1])
%     title(sprintf('Time step = %d',i));
%     xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
%     subplot(1,5,3)
%     mesh(x,y,gu3)
%     axis([-4 4 -4 4 0 1])
%     title(sprintf('Time step = %d',i));
%     xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
%     subplot(1,5,4)
%     mesh(x,y,gu4)
%     axis([-4 4 -4 4 0 1])
%     title(sprintf('Time step = %d',i));
%     xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
%     subplot(1,5,5)
%     mesh(x,y,gu5)
%     axis([-4 4 -4 4 0 1])
%     title(sprintf('Time step = %d',i));
%     xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
%     colormap cool
%     drawnow
% end

% curve graph
figure(2)
vu1=sum(sum(u1))/m^2; volu1(:)=vu1(1,1,:);
vu2=sum(sum(u2))/m^2; volu2(:)=vu2(1,1,:);
vu3=sum(sum(u3))/m^2; volu3(:)=vu3(1,1,:);
vu4=sum(sum(u4))/m^2; volu4(:)=vu4(1,1,:);
vu5=sum(sum(u5))/m^2; volu5(:)=vu5(1,1,:);
vu6=sum(sum(u6))/m^2; volu6(:)=vu6(1,1,:);
plot(t,volu1,'r',t,volu2,'b',t,volu3,'k',t,volu4,'g',t,volu5,'m',t,volu6,'LineWidth',2); grid
ylim([0 1])
legend('Prisoners','B or S','Chicken','Stag-Hunt','Hawk-Dove','generalized CK')
title('Replicator Dynamics in PDEs')
xlabel('time'); ylabel('cooperation rate');
