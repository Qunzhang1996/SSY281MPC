function Oinf
% Compute the maximum invariant set for an autonomous system using the MPT 
% toolbox

close all


%% Core code
% System Definition
% x(k+1) = A*x(k)
A=[0.5 0;1 -0.5];

% create model in MPT3 interface
model = LTISystem('A',A);

% constraints on states
model.x.min = [-10;-10];
model.x.max = [ 10; 10];

% constraint set
X = Polyhedron('lb',model.x.min,'ub',model.x.max);

% compute the maximum invariant set
S = model.invariantSet();


%% Plot results
% Printing parameters
label_font_size = 14;
tick_font_size  = 10;
line_width      = 0.8;
axeswidth       = 0.2;
set(0,'defaulttextinterpreter','latex');

figure
% plot the feasible set
plot(X,'color',[0.8 0.1 0.1],'linewidth',line_width);
hold on
% plot the invariant set
plot(S,'color',[0.2 0.4 0.6],'linewidth',line_width);
grid on

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
ht2=text(-1,-0.5,'$\mathcal{S}$');
set(ht2, 'FontSize', label_font_size);

border = 1;
axis([model.x.min(1)-border,model.x.max(1)+border,model.x.min(2)-border,model.x.max(2)+border])


end