function [induced_alpha,D_ind_2D,D_ind] = Induced_Drag_fun(config,Normals,L_2D,CPoint,InfiniteVortices, U_Inf, Gamma,Plot)

M = config.SemiSpanwiseDiscr;
N = config.ChordwiseDiscr;

induced_alpha = zeros(2*M(1),config.NBodies);       % induced angle of attack

% Compute induced alpha 
for kBody = 1:config.NBodies
    en = Normals{kBody}{1,1}.Coords;
    en = en / norm(en);
    for k = 1:2*M(kBody)
    V_vect = [0,0,0];
    control_point = CPoint{kBody,1}(k,:);
    for iBody = 1:config.NBodies
        for i = 1:2*M(iBody)
            for j = 1:N(iBody)
                % Gamma 
                gamma = Gamma{iBody,1}(j,i);

                % Compute the induced velocity by first  
                % semi-infinite vortex (at the root of the panel
                Extreme_a = InfiniteVortices{iBody,1}{j,i}.Root.toInfty;
                Extreme_b = InfiniteVortices{iBody,1}{j,i}.Root.onWing;
                U1 = induced_velocity_fun(control_point, Extreme_a, Extreme_b,gamma);

                % Compute the induced velocity by second            
                % semi-infinite vortex (at the tip of the panel)           
                Extreme_c = InfiniteVortices{iBody,1}{j,i}.Tip.toInfty;
                Extreme_d = InfiniteVortices{iBody,1}{j,i}.Tip.onWing;
                U2 = induced_velocity_fun(control_point, Extreme_d, Extreme_c,gamma);

                V_vect = V_vect + U1 + U2;
            end
        end
    end
    induced_alpha(k,kBody) = atan(dot(V_vect,en)/norm(U_Inf));
    end
end


% 2D Induced drag 
D_ind_2D = zeros(2*M(1),config.NBodies);
for iBody = 1:config.NBodies
    for i = 1:2*M(iBody)
    D_ind_2D(i,iBody) = abs(L_2D(i,iBody)*sin(induced_alpha(i,iBody)));
    end
end

% Total induced drag
D_ind_vect = zeros(2*M(1),config.NBodies);
for iBody = 1:config.NBodies
    for i = 1:2*M(iBody)
        delta_b = 2*config.SemiSpan(iBody)/(2*M(iBody));
        D_ind_vect(i,iBody) = D_ind_2D(i,iBody)*delta_b;
    end
end
if config.NBodies == 1
    D_ind = sum(D_ind_vect);
    fprintf('Total induced drag = %.3f N\n',D_ind);
elseif config.NBodies == 2
    D_ind_wing = sum(D_ind_vect(:,1));
    D_ind_tail = sum(D_ind_vect(:,2));
    D_ind = D_ind_wing + D_ind_tail;
    fprintf('Wing induced drag = %.3f N\n',D_ind_wing);
    fprintf('Tail induced drag = %.3f N\n',D_ind_tail);
    fprintf('Total induced drag = %.3f N\n',D_ind);
end

% Visualize induced angle of attack
if config.NBodies == 2 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    y2 = linspace(-config.SemiSpan(2),config.SemiSpan(2),2*M(2));
    plot(y1,rad2deg(induced_alpha(1:length(y1),1)),'o-','LineWidth',1.5);
    grid on 
    hold on 
    plot(y2,rad2deg(induced_alpha(1:length(y2),2)),'o-','LineWidth',1.5);
    title('Induced angles');
    xlabel('Spanwise position [m]');
    ylabel('[deg]');
    legend('Wing','Tail');
elseif config.NBodies == 1 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    plot(y1,rad2deg(induced_alpha(1:length(y1),1)),'o-','LineWidth',1.5);
    grid on 
    title('Induced angles');
    xlabel('Spanwise position [m]');
    ylabel('[deg]');
end

% Visualize induced drag per units of length
if config.NBodies == 2 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    y2 = linspace(-config.SemiSpan(2),config.SemiSpan(2),2*M(2));
    plot(y1,D_ind_2D(1:length(y1),1),'o-','LineWidth',1.5);
    grid on 
    hold on 
    plot(y2,D_ind_2D(1:length(y2),2),'o-','LineWidth',1.5);
    title('Induced Drag');
    xlabel('Spanwise position [m]');
    ylabel('[N/m]');
    legend('Wing','Tail');
elseif config.NBodies == 1 && Plot == 1
    figure
    y1 = linspace(-config.SemiSpan(1),config.SemiSpan(1),2*M(1));
    plot(y1,D_ind_2D(1:length(y1),1),'o-','LineWidth',1.5);
    grid on 
    title('Induced Drag');
    xlabel('Spanwise position [m]');
    ylabel('[N/m]');
end

end
