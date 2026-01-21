%% WEISSINGER METHOD 
close all 
clear all
clc

%% Test Case Sweep angle effect 
U_Inf_Mag = 62;     % [m/s]
alpha_values = [1, 2, 3];  % [deg]
beta = 0;           % [deg]
rho = 1.225;

config.NBodies = 1;

config.RootChord = [1];        % [m]
config.DihedralAngle = [0];    % [deg] 
config.SweepAngle = [0];       % [deg] (will be varied)
config.TaperRatio = [1];       
config.Span = [5];             % [m]

config.LEPosition_X = [0];     % [m]
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];     

config.RotationAngle_X = [0];       
config.RotationAngle_Y = [0];   
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [20];
config.ChordwiseDiscr = [20];

% Sweep angles to test
sweep_angles = [0, 10, 20]; % [deg]
results_CL = zeros(length(alpha_values), length(sweep_angles));
results_CDi = zeros(length(alpha_values), length(sweep_angles));

counter_alpha = 1;
for alpha = alpha_values
    counter_sweep = 1;
    for sweep_angle = sweep_angles
        config.SweepAngle = [sweep_angle];
        U_Inf = [cosd(alpha)*cosd(beta) sind(beta) sind(alpha)] .* U_Inf_Mag;

        %% Preliminary computations
        config.SemiSpan = config.Span./2;
        config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
        config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
        config.TipChord = config.RootChord .* config.TaperRatio;
        config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));
        M = config.SemiSpanwiseDiscr;
        N = config.ChordwiseDiscr;
        config.AspectRatio = config.Span.^2 ./ config.Surface;

        %% Create the geometry structure
        ControlPoints = cell(config.NBodies, 1);
        InducedPoints = cell(config.NBodies, 1);
        Normals = cell(config.NBodies, 1);
        InfiniteVortices = cell(config.NBodies, 1);
        Vortices = cell(config.NBodies, 1);
        internalMesh = cell(config.NBodies, 1); 
        WingExtremes = cell(config.NBodies, 1);

        for iBody = 1:config.NBodies
            [ControlPoints{iBody}, InducedPoints{iBody}, Normals{iBody}, InfiniteVortices{iBody}, Vortices{iBody}, internalMesh{iBody}, WingExtremes{iBody}] = createStructure(config, iBody);
        end

        %% Matrices initialization
        NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
        matrixA = zeros(NPanelsTot, NPanelsTot);
        knownTerm = zeros(NPanelsTot, 1);

        %% Construction of the matrix
        rowIndex = 0;
        for iBody = 1:config.NBodies
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    rowIndex = rowIndex + 1;
                    columnIndex = 0;
                    ControlPointHere = ControlPoints{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    for jCorpo = 1:config.NBodies
                        for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                            for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                                columnIndex = columnIndex + 1;
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                                U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                                Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                                Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                                matrixA(rowIndex, columnIndex) = dot(U, LocalNormal);
                            end
                        end
                    end
                end
            end
        end

        %% Known term construction 
        rowIndex = 0;
        for iBody = 1:config.NBodies
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    rowIndex = rowIndex + 1;
                    LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
                    knownTerm(rowIndex) = -dot(U_Inf, LocalNormal);
                end
            end
        end

        %% Solve the linear system
        Solution = linsolve(matrixA, knownTerm);
        Gamma = cell(config.NBodies, 1);
        rowIndex = 0;
        for iBody = 1:config.NBodies
            Gamma{iBody} = zeros( config.ChordwiseDiscr(iBody), config.SemiSpanwiseDiscr(iBody)*2 );
            for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
                    rowIndex = rowIndex + 1;
                    Gamma{iBody}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
                end
            end
        end

        %% Compute the 2D and 3D Lift
        [L_2D,L,delta_b] = Lift_fun(config,rho,U_Inf_Mag,Gamma,1);
        config.CL = L/(1/2 * rho * U_Inf_Mag^2 * config.Surface(1));
        results_CL(counter_alpha, counter_sweep) = config.CL;

        %% Compute 2D and 3D induced drag
        CPoint = C4_ControlPoints(config,internalMesh);
        [induced_alpha,D_ind_2D,D_ind] = Induced_Drag_fun(config,Normals,L_2D,CPoint,InfiniteVortices, U_Inf, Gamma,1);
        config.CD_ind = D_ind/(1/2 * rho * U_Inf_Mag^2 * config.Surface(1));
        results_CDi(counter_alpha, counter_sweep) = config.CD_ind;

        counter_sweep = counter_sweep + 1;
    end
    counter_alpha = counter_alpha + 1;
end

%% Display results
disp('CL and CDi for different sweep angles and angles of attack:');
disp(results_CL);
disp(results_CDi);

%% Plotting CL vs AoA for each sweep angle
figure;
hold on;
for i = 1:length(sweep_angles)
    plot(alpha_values, results_CL(:,i), 'DisplayName', ['Sweep = ' num2str(sweep_angles(i)) '°']);
end
grid on;
legend('Location', 'best');
xlabel('Angle of Attack (\alpha) [deg]');
ylabel('Lift Coefficient (C_L)');
title('C_L vs Angle of Attack for Different Sweep Angles');
hold off;

%% Optional: Plotting CDi vs AoA for each sweep angle
figure;
hold on;
for i = 1:length(sweep_angles)
    plot(alpha_values, results_CDi(:,i), 'DisplayName', ['Sweep = ' num2str(sweep_angles(i)) '°']);
end
grid on;
legend('Location', 'best');
xlabel('Angle of Attack (\alpha) [deg]');
ylabel('Induced Drag Coefficient (C_{Di})');
title('C_{Di} vs Angle of Attack for Different Sweep Angles');
hold off;

%% Visualization for last case
figure
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Aircraft Geometry')

colors = lines(config.NBodies);

for iBody = 1:config.NBodies
    for j = 1:config.ChordwiseDiscr(iBody)
        for i = 1:2*config.SemiSpanwiseDiscr(iBody)
            LE_root = internalMesh{iBody}{j,i}.LERoot;
            LE_tip  = internalMesh{iBody}{j,i}.LEtip;
            TE_tip  = internalMesh{iBody}{j,i}.TEtip;
            TE_root = internalMesh{iBody}{j,i}.TERoot;
            X = [LE_root(1) LE_tip(1) TE_tip(1) TE_root(1)];
            Y = [LE_root(2) LE_tip(2) TE_tip(2) TE_root(2)];
            Z = [LE_root(3) LE_tip(3) TE_tip(3) TE_root(3)];
            patch(X, Y, Z, colors(iBody,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k');
        end
    end
end

view(3)

% Color scale based on Gamma 
allGamma = [];
for iBody = 1:config.NBodies
    allGamma = [allGamma; Gamma{iBody}(:)];
end
cmin = min(allGamma);
cmax = max(allGamma);
colormap(jet)
clim([cmin cmax])
cb = colorbar;
cb.Label.String = '\Gamma [m^2/s]';

% Plot bound vortices
for iBody = 1:config.NBodies
    for j = 1:config.ChordwiseDiscr(iBody)
        for i = 1:2*config.SemiSpanwiseDiscr(iBody)
            V = Vortices{iBody}{j,i};
            x = [V.Root(1), V.Tip(1)];
            y = [V.Root(2), V.Tip(2)];
            z = [V.Root(3), V.Tip(3)];
            gamma_ji = Gamma{iBody}(j,i);
            t = (gamma_ji - cmin) / (cmax - cmin);
            t = max(0,min(1,t));
            color = interp1(linspace(0,1,64), jet(64), t);
            plot3(x,y,z,'LineWidth',2,'Color',color);
        end
    end
end

save config 
view(3)
