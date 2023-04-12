%% Main script of PSS1 in SSY281
close all

% Our test system
A  = [1.0025 0.1001;0.05 1.0025];
B  = [0.005;0.1001];
C  = [1 0];
D  = 0;
Q  = [5 0;0 1];
Pf = [5 0;0 1];
R  = 0.5;
n  = length(A); % 2
m  = size(B,2); % 1

N  = 10;    % set horizon length               
x0 = [1;0]; % initial state


%% Question 1 - Unconstrained RHC (LQ problem)
[Z,VN] = URHC(A,B,N,Q,R,Pf,x0,n); % solve QP numerically 
x      = Z(1:n*N);                % get optimal state sequence
u      = Z(n*N+1:end);            % get optimal control input sequence
u0     = u(1);                    % control applied by the RHC

pause


%% Question 2 - Constrained RHC
% Assemble constraints in matrix form as: Fx + Gu <= h 
% Here, assume the box constraints: abs(x2(k)) <= x2_max and abs(u(k)) <= u_max for all k  
x2_max = 0.5; % start with 0.5 (but feel free to play around with this)
u_max  = 0.7; % start with 0.7 (but feel free to play around with this)

% Fill this part.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F      = kron([eye(N); zeros(N)], [0 1;0 -1]);
G      = kron([zeros(N); eye(N)], [1;-1]);
h      = [x2_max*ones(n*N,1); u_max*ones(n*N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Z,VN] = CRHC(A,B,N,Q,R,Pf,F,G,h,x0,n); % solve QP numerically
x      = Z(1:n*N);                      % get optimal state sequence
u      = Z(n*N+1:end);                  % get optimal control input sequence
u0     = u(1);                          % control applied by the RHC

pause


%% Question 3 - Time-domain simulations
% Change this option to change which type of plot you want
plot_option = 2; % 1 for time domain plots and 2 for x1-x2 space plots

% Simulation parameters (feel free to change tf)
tf     = 50;               % assume a sampling interval of 1
x_vec  = [x0 zeros(n,tf)]; % store state trajectory
u_vec  = zeros(m,tf);      % store control input trajectory
x      = x0;               % initial condition

% Carry out simulations by applying the receding horizon principle
% Fill this part.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:tf
    % Solve RHC problem (either unconstrained or constrained)
    [Z,~] = CRHC(A,B,N,Q,R,Pf,F,G,h,x,n);
    
    % Get control input to apply and next state
    x = Z(1:n);
    u = Z(n*N+1);

    % Save to structure
    x_vec(:,iter+1) = x;
    u_vec(:,iter)   = u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot figures (you don't need to worry about this) 
% Time domain plots 
if plot_option == 1
    % Plot state and command trajectories
    clf;
    subplot(3,1,1); hold on;
    plot(0:tf,x_vec(1,:),'Linewidth',2)
    title('State x1')

    subplot(3,1,2); hold on;
    plot(0:tf,x_vec(2,:),'r','Linewidth',2)
    plot(0:tf,-x2_max*ones(1,tf+1),'k--')
    plot(0:tf,x2_max*ones(1,tf+1),'k--')
    title('State x2')

    subplot(3,1,3); hold on;
    plot(0:tf-1,u_vec,'m','Linewidth',2)
    plot(0:tf-1,-u_max*ones(1,tf),'k--')
    plot(0:tf-1,u_max*ones(1,tf),'k--')
    title('Control input u')
% x1-x2 space plots
elseif plot_option == 2
    % Plot state trajectory progressively 
    for iter = 1:tf+1
        clf; hold on;
        scatter(x_vec(1,1:iter),x_vec(2,1:iter),'filled')      % trajectory
        scatter(0,0,'rx')                                      % origin
        plot(linspace(-0.2,1.5,100),-x2_max*ones(1,100),'k--') % upper bound 
        plot(linspace(-0.2,1.5,100),x2_max*ones(1,100),'k--')  % lower bound
        xlim([-0.2 1.5]) % fix x and y range to avoid constant rescaling
        ylim([-x2_max-0.1 x2_max+0.1])
        pause(0.05);     % briefly pause plotting process (feel free to comment this line out)
    end
else 
    disp('Plot option not recognized')
end



