%% ODE parameter : HD
clear, clc
subplot(1,3,2)
b=0:0.01:3; V=b; C=1;
for i=1:length(V)
    Payoff_HD=[(V(i)-C)/2 V(i); 0 V(i)/2];
    dxdt =@(t,x) [x(1)*((V(i)-C)/2*x(1)*(1-x(1))+V(i)*x(2)*(1-x(1)-x(2)/2)); x(2)*(V(i)/2*x(2)-(V(i)-C)/2*x(1)^2-V(i)*x(1)*x(2)-V(i)/2*x(2)^2)];
    [t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
    Y(i)=sum(y(:,1))/length(y(:,1));
end
plot(V,Y,'mo');grid
xlabel('Temptation'); ylabel('Mean of the cooperators')
ylim([-0.1 1.1])
hold on
%% PDE parameter : HD
clear, clc
rng('default')
% parameters
b=0:0.01:3; c=1; v=b; fr=(v-c)/2; h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C=zeros(m,l,n);D=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C(:,:,1)=ones(m,l)*0.9+r*0.4;
D(:,:,1)=ones(m,l)*0.1+r*0.4;
for z=1:length(v)
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
                    dt*C(j,k,i)*(fr(z)*C(j,k,i)*(1-C(j,k,i))+v(z)*D(j,k,i)*(1-C(j,k,i)-D(j,k,i)/2));
                % Defectors
                D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                    dt*D(j,k,i)*(v(z)/2*D(j,k,i)-fr(z)*C(j,k,i)^2-v(z)*C(j,k,i)*D(j,k,i)-v(z)/2*D(j,k,i)^2);
            end
        end
    end
    vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
    Y(z)=sum(volC)/length(volC);
end
plot(b,Y,'b-o');
%% ------------------------------------------------------------
clear
clc
% Index
Col = [1 2
    3 4];
% Each row defines: red, green, yellow and blue
map = [ 0 0 1
    0 1 0
    1 1 0
    1 0 0];
fid = 1; %Output directed to screen (=1) or elsewhere...
sep ='-----------';
% Maximum number of iterations
max_iter = 100;
% Checkboard size
N = 100;
% Number of cells
n = N*N;
% Proportion of defectors
p = 0.1;
% Payoff matrix
% Payoff = [-0.3 1.9; 0 0.95] % HD
% Payoff = [0.45 -0.05; 0.95 0] % SD
B=0:0.01:3; V=B; C=1;
for z=1:length(B)
    Payoff=[(V(z)-C)/2 V(z); 0 V(z)/2]; % HD
    % Allocations
    pC = zeros(max_iter,1);
    Payment = zeros(N,N);
    NE = zeros(N,N);
    Tab = zeros(N,N);
    Tab_cube = zeros(N,N,max_iter);
    % ------------------------------------------------------------
    % INITIAL CONDITIONS
    E = ones(N,N);
    % % Random
    A = rand(1,N*N);
    % % Changing E according to the initial proportion of defectors
    I = find(A < p);
    E(I) = 2;
    
    Tab_cube(:,:,1)=E(:,:);
    pC(1)=length(find(E==1))/n; pD(1)=length(find(E==2))/n;
    
    for iter=2:max_iter
        % Setting payments when each player plays with their 8 neighbours,
        % including herself.
        for i=1:N
            for j=1:N
                pa = 0;
                for k=-1:1
                    for h=-1:1
                        % Taking account of boundary conditions
                        a = cdc(i+k,N); b = cdc(j+h,N);
                        % Setting payment accordint to the selected strategy
                        pa = pa + Payoff(E(i,j), E(a,b));
                    end
                end
                Payment(i,j) = pa;
            end
        end
        % Evaluation of the environment and possible change of strategy
        for i=1:N
            for j=1:N
                pay = Payment(i,j);
                NE(i,j) = E(i,j);
                for k=-1:1
                    for h=-1:1
                        % Taking account of boundary conditions
                        a = cdc(i+k,N); b = cdc(j+h,N);
                        % If the neighbour performed better, the i,j player
                        % changes and mimics the neighbour
                        if (Payment(a,b) > pay)
                            pay = Payment(a,b);
                            NE(i,j) = E(a,b);
                        end
                    end
                end
            end
        end
        %
        % Changes of strategy becomes effective here.
        for i=1:N
            for j=1:N
                Tab(i,j) = Col(NE(i,j),E(i,j));
                % Change of strategy (updating matrix E)
                E(i,j) = NE(i,j);
            end
        end
        Tab_cube(:,:,iter) = Tab;
        % Proportion of ccoperators
        pC(iter) = length(find(Tab==1 | Tab==2)) / n;
        pD(iter) = length(find(Tab==3 | Tab==4)) / n;
    end %iter
    Y(z)=sum(pC)/length(pC);
end
plot(V,Y,'g-o'); grid
hold off