clc; clear all;
set(0,'defaulttextinterpreter','latex');
label_font_size = 14;


%% Defining system
model = LTISystem('A', [1.5 0;1 -1.5], 'B', [1; 0]);     % unstable A
% model = LTISystem('A', [0.5 0.9; 0 0.9], 'B', [0; 1]); % stable A
model.x.min = [-10; -10];
model.x.max = [ 10;  10];
model.u.min = -5;
model.u.max =  5;

X = Polyhedron('lb', model.x.min, 'ub', model.x.max);
U = Polyhedron('lb', model.u.min, 'ub', model.u.max);


%% Defining Controller Tuning    
Q = [1 0; 0 1];
model.x.penalty = QuadFunction(Q);
R = 1;
model.u.penalty = QuadFunction(R);
N = 5;

P      = model.LQRPenalty; % i.e. P is the DARE solution for this LQR
LQRset = model.LQRSet;     % largest control invariant set for closed-loop LQR system
C_inf  = model.invariantSet(); % largest control invariant set

% Several possible choices of terminal set
% Tset   = C_inf;
% Tset   = LQRset;
Tset   = Polyhedron([0 0]);

model.x.with('terminalPenalty');
model.x.terminalPenalty = P;

model.x.with('terminalSet');
model.x.terminalSet = Tset;

% MPC and EMPC
tic
mpc  = MPCController(model, N);
toc
tic
empc = mpc.toExplicit();
toc

% Compute the feasible set XN
XN = Tset;
for i = 1:N
    XN = model.reachableSet('X', XN, 'U', U, 'direction', 'backward');
    XN = XN.intersect(X).minHRep();
end


%% Plot sets
figure(1)
plot(X, 'color', 'r', C_inf, 'color', 'b', XN, 'color', 'g', Tset, 'color', 'y');

% Labels
hx1 = xlabel('$x_1$');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$');
set(hy1, 'FontSize', label_font_size);

% Legend
ht1=text(9,9,'$\mathbf{X}$');
set(ht1, 'FontSize', label_font_size);
ht2=text(2,-3,'$\mathcal{C}_{\infty}$');
set(ht2, 'FontSize', label_font_size);
ht3=text(-2,3,'$\mathcal{X}_N$');
set(ht3, 'FontSize', label_font_size);
ht4=text(0.5,0.5,'$\mathbf{X}_f$');
set(ht4, 'FontSize', label_font_size);

% Plot EMPC critical regions 
figure(2)
empc.partition.plot(); % Plot the explicit MPC partition (2D)
% empc.feedback.fplot();   % Plot the partition in 3D (u is on the z-axis)

hx1 = xlabel('$x_1$');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$');
set(hy1, 'FontSize', label_font_size);


%% Compute and plot closed-loop simulations
Nsim = 15; % number of iterations
loop  = ClosedLoop(mpc, model);
eloop = ClosedLoop(empc, model);

% Choose the starting point (e.g. one vertex of XN)
x0 = XN.V(1,:)';
% x0 = [-2.3;-3.6]; % alternative starting point

% Simulate (compare time difference)
tic
data  = loop.simulate(x0, Nsim);
toc
tic
edata = eloop.simulate(x0, Nsim);
toc


% Plot state trajectories
figure(3); 
clf; 
subplot(2,1,1); 
hold on;
plot(0:Nsim,data.X(1,:), 'Linewidth', 2)
plot(0:Nsim,data.X(2,:), 'Linewidth', 2)
stairs(data.U, 'k--')
xlabel('Iteration')
title('MPC')
legend({'$x_1$','$x_2$','$u$'}, 'interpreter', 'latex', 'location', 'southeast')

subplot(2,1,2); 
hold on;
plot(0:Nsim,edata.X(1,:), 'Linewidth', 2)
plot(0:Nsim,edata.X(2,:), 'Linewidth', 2)
stairs(edata.U, 'k--')
xlabel('Iteration')
title('EMPC')
legend({'$x_1$','$x_2$','$u$'}, 'interpreter', 'latex')
