clear
close all
clc


H = readmatrix ('cp_bubbledancer.txt');
CpX = H(:,2);
X = H(:,1);
plot(X, CpX, 'r', 'DisplayName', 'Text file data');
xlabel('x/c');
ylabel('C_p');
title('C_p distribution');
set(gca,'YDir','reverse'); % We print -Cp
grid on;
hold on;
H = readmatrix ('cp0008_re500000.txt');
CpX = H(:,2);
X = H(:,1);
plot(X, CpX, 'b', 'DisplayName', 'Text file data');
xlabel('x/c');
ylabel('C_p');
title('C_p distribution');
set(gca,'YDir','reverse'); % We print -Cp
legend('AG26 Bubble Dancer','NACA0008','Text file data','Location','best');
grid on;


% Transition point bubble dancer: 0.4887(x/c)
% Transition point NACA 0008: 0.3521(x/c)