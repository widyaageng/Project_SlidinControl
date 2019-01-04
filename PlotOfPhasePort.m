%%--------------------------
%----- SLIDING CONTROL------
%%--------------------------

%% Plot of Phase Portrait

x1 = state_vars.signals.values(:,1);
x2 = state_vars.signals.values(:,2);
time = state_vars.time;
surf_line  = surf_C(2)*(-1)*time;

figure;
hold on;
plot(time,x1,'-.r');
plot(time,x2,'-.b');
plot(time,surf_line,'-g');
legend('x1','x2','surf');
axis([time(1) time(end)+1 min([min(x1) min(x2) min(surf_line)])-0.5 max([max(x1) max(x2) max(surf_line)])+0.5]);
grid; hold off;