function PreReachAut
% Compute the one-step forward and backward reachable sets using the MPT
% toolbox for autonomous systems

close all


%% Core code
% Autonomous system
% x(k+1) = A*x(k)
A=[0.5 0;1 -0.5]; % asymptotically stable

% create model in MPT3 interface
model = LTISystem('A',A);

% constraints on states
model.x.min = [-10;-10];
model.x.max = [ 10; 10];

% constraint sets represented as polyhedra
X = Polyhedron('lb',model.x.min,'ub',model.x.max);


%% Backward reachable set (preset)
Xpre = model.reachableSet('X', X, 'direction', 'backward','N',1);

% intersect with the state constraints
XpreandX = Xpre.intersect(X).minHRep();


%% Forward reachable set
Xreach1 = model.reachableSet('X', X, 'direction', 'forward','N',1);


%% Plot results
% Printing parameters
label_font_size = 14;
tick_font_size  = 10;
line_width      = 0.8;
axeswidth       = 0.2;
set(0,'defaulttextinterpreter','latex');

fig1 = figure;
% plot the set X
plot(X,'color',[0.8 0.1 0.1],'linewidth',line_width);
hold on
% plot the preset of X
% plot(Xpre,'color',[0.2 0.4 0.6],'linewidth',line_width); 
plot(XpreandX,'color',[0.2 0.4 0.6],'linewidth',line_width); 
border = 1;
axis([model.x.min(1)-border,model.x.max(1)+border,model.x.min(2)-border,model.x.max(2)+border])

ht1=text(7,-8,'$\mathcal{X}$');
set(ht1, 'FontSize', label_font_size);
ht2=text(-6,-0.5,'Pre$(\mathcal{X})\cap\mathcal{X}$');
set(ht2, 'FontSize', label_font_size);
xt = transpose(-10:5:10);
yt = transpose(-10:5:10);
set(gca,'XTick',xt);
set(gca,'YTick',yt);
set(gca,'YTickLabel',num2str(xt));
set(gca,'XTickLabel',num2str(yt));

set(gca,'LineWidth',axeswidth)
set(gca,'FontSize', tick_font_size);

hx1 = xlabel('$x_1$');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$');
set(hy1, 'FontSize', label_font_size);


fig2 = figure;
% plot the set X
plot(X,'color',[0.8 0.1 0.1],'linewidth',line_width);
hold on
% plot forward reachable set
plot(Xreach1,'color',[0.2 0.4 0.6],'linewidth',line_width); 
border2 = 6;
axis([model.x.min(1)-border2,model.x.max(1)+border2,model.x.min(2)-border2,model.x.max(2)+border2]);

set(gca,'LineWidth',axeswidth)
set(gca,'FontSize', tick_font_size);

hx1 = xlabel('$x_1$');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$');
set(hy1, 'FontSize', label_font_size);

ht3=text(7,-8,'$\mathcal{X}$');
set(ht3, 'FontSize', label_font_size);
ht4=text(-3.5,-7,'Suc($\mathcal{X}$)');
set(ht4, 'FontSize', label_font_size, 'Rotation', 57);

xt = transpose(-10:5:10);
yt = transpose(-10:5:10);
set(gca,'XTick',xt);
set(gca,'YTick',yt);
set(gca,'YTickLabel',num2str(yt));
set(gca,'XTickLabel',num2str(xt));


end
