%% ODE parameter : PD
clear, clc
b=0:0.01:3;
subplot(1,3,1)
% Y0=zeros(length(b));
% Y=Y0(1,:);
for i=1:length(b)
    Payoff_PD=[1 0; b(i) 0];
    dxdt =@(t,x) [x(1)^2*(1-x(1)-b(i)*x(2)); x(1)*x(2)*(b(i)*(1-x(2))-x(1))];
    [t,y]=ode45(dxdt,[0 100],[0.9 0.1]);
    Y(i)=sum(y(:,1))/length(y(:,1));
end
plot(b,Y,'ro');grid
xlabel('Temptation'); ylabel('Mean of the cooperators')
ylim([-0.1 1.1])
hold on

%% PDE parameter : PD
clear, clc
rng('default')
% parameters
b=0:0.01:3; h=2/sqrt(10); dt=(h^2)/5;
t=0:dt:100; x=0:h:100; y=0:h:100;
n=length(t); m=length(x); l=length(y);
C=zeros(m,l,n);D=zeros(m,l,n);

% initial conditions
r = -1 + (1+1)*rand(m,l);
C(:,:,1)=ones(m,l)*0.9+r*0.4;
D(:,:,1)=ones(m,l)*0.1+r*0.4;

for z=1:length(b)
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
                    dt*C(j,k,i)^2*((1-C(j,k,i))-b(z)*D(j,k,i));
                % Defectors
                D(j,k,i+1)=D(j,k,i)+dt*((D(j+1,k,i)+D(j-1,k,i)+D(j,k+1,i)+D(j,k-1,i)-4*D(j,k,i))/h^2)+...
                    dt*C(j,k,i)*D(j,k,i)*(b(z)*(1-D(j,k,i))-C(j,k,i));
            end
        end
    end
    
    % graphs
    vC=sum(sum(C))/m^2; volC(:)=vC(1,1,:);
    % vD=sum(sum(D))/m^2; volD(:)=vD(1,1,:);
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
B=0:0.01:3;
for z=1:length(B)
    Payoff=[1 0; B(z) 0]; % PD
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
plot(B,Y,'g-o'); grid
hold off

