function [L_2D,L,delta_b] = Lift_fun(config,rho,U_Inf_Mag,Gamma,Plot)
% This functions computes the 2D lift of each section and the total lift of
% a 3D wing for the Weissinger method
% NOTE: Plot == 1 means that there will be plots
%       Plot != 1 means no plots

% 2D lift - Single Wing
M = config.SemiSpanwiseDiscr;
L_2D = zeros(2*M(1),config.NBodies);
for iBody = 1:config.NBodies
    for i = 1:2*M(iBody)
    L_2D(i,iBody) = sum(rho*U_Inf_Mag*Gamma{iBody}(:,i)*cosd(config.DihedralAngle(iBody)));
    % Gamma{body}(chord position,span position)
    end
end

% 3D Lift - Single wing
L_vect = zeros(2*M(1),config.NBodies);
for iBody = 1:config.NBodies
    for i = 1:2*M(iBody)
        delta_b = 2*config.SemiSpan(iBody)/(2*M(iBody));
        L_vect(i,iBody) = L_2D(i,iBody) * delta_b;
    end
end

% Values of lift
if config.NBodies == 1
    L_wing = sum(L_vect(:,1));
    L = sum(sum(L_vect));
    fprintf('Wing Lift = %.3f N\n',L_wing);
    fprintf('Total Lift = %.3f N\n',L);
elseif config.NBodies == 2
    L_wing = sum(L_vect(:,1));
    L_tail = sum(L_vect(1:2*M(2),2));
    L = sum(sum(L_vect));
    fprintf('Wing Lift = %.3f N\n',L_wing);
    fprintf('Tail Lift = %.3f N\n',L_tail);
    fprintf('Total Lift = %.3f N\n',L);
end

if config.NBodies == 2 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    y2 = linspace(-config.SemiSpan(2),config.SemiSpan(2),2*M(2));
    plot(y1,L_2D(1:length(y1),1),'o-','LineWidth',1.5);
    grid on 
    hold on 
    plot(y2,L_2D(1:length(y2),2),'o-','LineWidth',1.5);
    title('2D Lift distribution - Wing and tail');
    xlabel('Spanwise position [m]');
    ylabel('[N/m]')
    legend('Wing','Tail');
elseif config.NBodies == 1 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    plot(y1,L_2D(1:length(y1),1),'o-','LineWidth',1.5);
    grid on 
    title('2D Lift distribution - Wing and tail');
    xlabel('Spanwise position [m]');
    ylabel('[N/m]')
end
end

