function Cinf
% Compute the maximum control invariant set using the MPT toolbox

close all


%% Core code
% matrices of LTI dynamics 
% x(k+1) = A*x(k) + B*u(k)
A=[1.5 0;1 -1.5];
B=[1;0];

% create model in MPT3 interface
model = LTISystem('A',A,'B',B);

% constraints on inputs and states
model.u.min = -5;
model.u.max = 5;
model.x.min = [-10;-10];
model.x.max = [ 10; 10];

% constraint sets represented as polyhedra
X = Polyhedron('lb',model.x.min,'ub',model.x.max);
U = Polyhedron('lb',model.u.min,'ub',model.u.max);

% compute the invariant set including the intermediate steps
maxIterations = 10;
Piter         = []; % save intermediate set during computation
P            = X;  % initial set constraint
% recursive computation of the maximum invariant set
for i = 1:maxIterations
    %%%%%%% Fill this part %%%%%%%%%%%
    % It should return the maximum
    % control invariant set P, and the
    % sequence of intermediate sets
    % in Piter.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Omega(i+1) = pre(Omega(i) and Omega(i)
    Piter = [Piter P];
    back_P = model.reachableSet('X', P, 'U', U, 'direction', 'backward','N',1);
    P = P.intersect(back_P).minHRep();

end

% direct command for computing the invariant set (but here we want to show
% intermediate sets)
% Cinf = model.invariantSet(); 

% the invariant set
Cinf = P;


%% Plot results
% Printing parameters
label_font_size = 14;
tick_font_size  = 10;
line_width      = 0.8;
axeswidth       = 0.2;
set(0,'defaulttextinterpreter','latex');

figure
% plot the constraint set
plot(X,'color',[0.8 0.1 0.1],'linewidth',line_width);
hold on
grid on
%% Students, here you can plot the intermediate sets
for i = 1:length(Piter)
    plot(Piter(i),'color',[0.2 0.4 0.6],'alpha',i/length(Piter))
end
%%
% plot the invariant set
plot(Cinf,'color',[0.2 0.4 0.6],'linewidth',line_width)
axis([model.x.min(1),model.x.max(1),model.x.min(2),model.x.max(2)])

set(gca,'LineWidth',axeswidth)
set(gca,'FontSize', tick_font_size);
xt = transpose(-10:5:10);
yt = transpose(-10:5:10);
set(gca,'XTick',xt);
set(gca,'YTick',yt);
set(gca,'YTickLabel',num2str(xt));
set(gca,'XTickLabel',num2str(yt));

hx1 = xlabel('$x_1$');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$');
set(hy1, 'FontSize', label_font_size);

ht1=text(7,-8,'$\mathcal{X}$');
set(ht1, 'FontSize', label_font_size);
ht2=text(-1,-0.5,'$\mathcal{C}_{\infty}$');
set(ht2, 'FontSize', label_font_size);

border = 1;
axis([model.x.min(1)-border,model.x.max(1)+border,model.x.min(2)-border,model.x.max(2)+border])


end