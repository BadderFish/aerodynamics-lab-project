%% Comparison between Weissinger and XFLR5 
clear 
close all
clc

%% Data loading and error estimations: SINGLE WING
% Load Weissinger Data
load Wessinger.mat
alpha_Weis = 0:1:5;
CL_Weis = config.CL;
CDi_Weis = config.CD_ind;

% Load XFLR5 Data
CL_alpha = load('CL-alpha.dat');
CDi_alpha = load('CDi-alpha.dat');
alpha_XFLR5 = CL_alpha(:,1);
CL_XFLR5 = CL_alpha(:,2);
CDi_XFLR5 = CDi_alpha(:,2);

% Error calculation 
CL_error = abs(CL_XFLR5 - CL_Weis)./CL_XFLR5 * 100;
CDi_error = abs(CDi_XFLR5 - CDi_Weis)./CDi_XFLR5 * 100;

%% Plots 
% CL plot
figure
sgtitle('Weissinger - XFLR5 Comparison: Single Wing');
subplot(2,2,1)
plot(alpha_Weis,CL_Weis,'o-','LineWidth',1.5);
grid on 
hold on 
plot(alpha_XFLR5,CL_XFLR5,'o-','LineWidth',1.5);
title('CL_{\alpha} curve');
xlabel('\alpha [deg]');
ylabel('CL');
legend('Weissinger','XFLR5');

subplot(2,2,2);
plot(alpha_XFLR5,CL_error,'o-','LineWidth',1.5);
grid on 
title('Relative error of CL');
xlabel('\alpha [deg]');
ylabel('Relative error [%]');

% CDi plots

sgtitle('Weissinger - XFLR5 Comparison: Single Wing');
subplot(2,2,3)
plot(alpha_Weis,CDi_Weis,'o-','LineWidth',1.5);
grid on 
hold on 
plot(alpha_XFLR5,CDi_XFLR5,'o-','LineWidth',1.5);
title('CDi_{\alpha} curve');
xlabel('\alpha [deg]');
ylabel('CDi');
legend('Weissinger','XFLR5');

subplot(2,2,4);
plot(alpha_XFLR5,CDi_error,'o-','LineWidth',1.5);
grid on 
title('Relative error of CDi');
xlabel('\alpha [deg]');
ylabel('Relative error [%]');

clearvars 

%% Data loading and error estimations: COMPLETE AIRCRAFT
% Load Weissinger Data
load Wessinger_complete.mat
alpha_Weis = 1:1:5;
CL_Weis = config.CL(2:end);
CDi_Weis = config.CD_ind(2:end);
% Data from XFLR5 only from 1 to 5 degrees 

% Load XFLR5 Data
CL_alpha = load('CL-alpha_complete.dat');
CDi_alpha = load('CDi-alpha_complete.dat');
alpha_XFLR5 = CL_alpha(:,1);
CL_XFLR5 = CL_alpha(:,2);
CDi_XFLR5 = CDi_alpha(:,2);

% Error calculation 
CL_error = abs(CL_XFLR5 - CL_Weis)./CL_XFLR5 * 100;
CDi_error = abs(CDi_XFLR5 - CDi_Weis)./CDi_XFLR5 * 100;

%% Plots 
% CL plot
figure
sgtitle('Weissinger - XFLR5 Comparison: Wing + Horizontal Tail');
subplot(1,2,1)
plot(alpha_Weis,CL_Weis,'o-','LineWidth',1.5);
grid on 
hold on 
plot(alpha_XFLR5,CL_XFLR5,'o-','LineWidth',1.5);
title('CL_{\alpha} curve');
xlabel('\alpha [deg]');
ylabel('CL');
legend('Weissinger','XFLR5');

subplot(1,2,2);
plot(alpha_XFLR5,CL_error,'o-','LineWidth',1.5);
grid on 
title('Relative error of CL');
xlabel('\alpha [deg]');
ylabel('Relative error [%]');

% CDi plots
figure
sgtitle('Weissinger - XFLR5 Comparison: Wing + Horizontal Tail');
subplot(1,2,1)
plot(alpha_Weis,CDi_Weis,'o-','LineWidth',1.5);
grid on 
hold on 
plot(alpha_XFLR5,CDi_XFLR5,'o-','LineWidth',1.5);
title('CDi_{\alpha} curve');
xlabel('\alpha [deg]');
ylabel('CDi');
legend('Weissinger','XFLR5');

subplot(1,2,2);
plot(alpha_XFLR5,CDi_error,'o-','LineWidth',1.5);
grid on 
title('Relative error of CDi');
xlabel('\alpha [deg]');
ylabel('Relative error [%]');