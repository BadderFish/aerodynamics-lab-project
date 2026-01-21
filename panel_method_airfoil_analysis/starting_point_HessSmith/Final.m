%% test Hess Smith

clc
close all
clear 

addpath mat_functions
addpath Xfoil_cp_distributions
%% Input

U_inf = 1;      % Far-field velocity [m/s]
AoA = -4;       % Angle of attack (you have to take it negative because it is the free stream angle of attack)
U_inf_x = U_inf * cos(deg2rad(AoA));
U_inf_y = U_inf * sin(deg2rad(AoA));
U_inf = [U_inf_x; U_inf_y];

Chord = 1;
NPanels = 301;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Create profile (with xfoil)

[x,y]=createProfile('0008_301',NPanels,Chord);

% You just need to name it like this, changing the number of panels, and
% you are done

geo.x=x;
geo.y=y;

figure
plot(x,y,'o-')
axis equal


%% Create discretization & initialization

[centers,normals,tangent,extrema_1,extrema_2,alpha,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = CreatePanels(geo);
        
NCols = sum(NPanels) + 1;
NRows = NCols;
A = zeros(NRows,NCols);     % system coefficients
B = zeros(NRows,1);         % known terms

%% Fill A

for i = 1:NPanels
    local_center = centers(i, :)';
    local_normal = normals(i, :)';

    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';

        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

        A(i, j) = dot(uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);

        A(i, sum(NPanels)+1) = A(i, sum(NPanels)+1) + dot(uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);
    end

end

%% Create a_v, c_s, and c_v vectors

first_centers = centers(1, :)';
first_tangent = tangent(1, :)';

last_centers = centers(end, :)';
last_tangent = tangent(end, :)';

last_a = 0;
for j = 1:NPanels
    local_extreme_1 = extrema_1(j, :)';
    local_extreme_2 = extrema_2(j, :)';
    local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
    local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(uSource(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);
    last_a = last_a + dot(uVortex(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);

    a = a + dot(uSource(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);
    last_a = last_a + dot(uVortex(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);

    A(sum(NPanels) + 1, j) = a;
end

A(sum(NPanels) + 1, sum(NPanels) + 1) = last_a;

%% Create B, the known terms in the system 

for j = 1:NPanels
    local_normal = normals(j, :)';
    B(j) = - dot(U_inf, local_normal);
end

first_tangent = tangent(1, :)';
last_tangent = tangent(end, :)';
B(sum(NPanels) + 1) = - dot(U_inf, (first_tangent + last_tangent));

%% Solve the linear system

solution = linsolve(A,B);

%% Compute velocity on control points and pressure distribution

q     = solution(1:NPanels);
gamma = solution(NPanels+1);

vel_vec = zeros(NPanels,2);   % velocity vectors (x,y) in control points
u_t     = zeros(NPanels,1);   % tangential component
Cp      = zeros(NPanels,1);   % pressure coefficient

Uinf_vec = U_inf;
Uinf_mag = norm(Uinf_vec);
Uinf_unit = Uinf_vec / Uinf_mag;

% Compute velocity in each i-control point
for i = 1:NPanels
    xc = centers(i, :)';        % control point (2x1)
    ui = Uinf_vec;              % 

    for j = 1:NPanels
        ex1 = extrema_1(j, :)'; % extrema of panel j
        ex2 = extrema_2(j, :)';
        L2G = squeeze(L2G_TransfMatrix(j, :, :));
        G2L = squeeze(G2L_TransfMatrix(j, :, :));

        sou = uSource(xc, ex1, ex2, L2G, G2L);   
        vor = uVortex(xc, ex1, ex2, L2G, G2L);   

        ui = ui + q(j)*sou + gamma*vor;
    end

    vel_vec(i, :) = ui';                        % store (x,y)
    
    % Compute tangential component
    ti = tangent(i, :)';                        
    u_t(i) = dot(ui, ti);
    Cp(i) = 1 - (u_t(i)/Uinf_mag)^2;            % Incompressible Bernoulli
end

%% Circulation, Cl (Kutta-Joukowsky) and Cl integrating pressure

Gamma_total = gamma * sum(lengths);    

% Unit vector perpendicular to freestream (rotate CCW by 90deg)
n_perp = [-Uinf_unit(2); Uinf_unit(1)];       

% Projection of panel normals onto lift direction (perp to U_inf)
n_perp_proj = normals * n_perp;   % Nx1 vector: n · n_perp

Cl_pressure = sum( Cp .* lengths' ./ Chord .* n_perp_proj );

% Kutta-Joukowsky (dimensionless) (for comparison)
Cl_KJ = -2 * Gamma_total / (Uinf_mag * Chord);

% Dimensional: L = rho * U * Gamma_total
rho = 1.225;   % density at 0 height
L_KJ = -rho * Uinf_mag * Gamma_total;
Cl_from_L = L_KJ / (0.5 * rho * Uinf_mag^2 * Chord); 

%% Aerodynamic moment wrt Leading Edge (CM,LE)
%LE_coord = zeros(NPanels,2);
LE_coord = repmat([LE_X_Position, LE_Y_Position], NPanels, 1);
rc = centers - LE_coord;

cross_z = rc(:,1).*normals(:,2) - rc(:,2).*normals(:,1);
CM_LE = - sum( Cp .* lengths'./Chord .* cross_z );

CM_AC = CM_LE + 0.25*Cl_KJ;

%% Results
fprintf('U_inf = %.6f m/s\n', Uinf_mag);
fprintf('alpha = %.6f °\n', AoA);
fprintf('Gamma total = %.6e\n', Gamma_total);
fprintf('C_l (pressure integral) = %.6f\n', Cl_pressure);
fprintf('C_l (Kutta-Joukowsky)  = %.6f\n', Cl_KJ);
fprintf('C_m,LE = %.6f\n', CM_LE);
fprintf('C_m,AC = %.6f\n', CM_AC);




%% Plot

x_over_c = centers(:,1) / Chord;

% 1. Find the index of the leading edge (min x)
[~, idx_LE] = min(x_over_c);

% 2. Separate into Lower and Upper Surface indices
% Assuming panels go: Lower (1 to idx_LE) -> Upper (idx_LE to NPanels)

% Lower Surface (usually from TE to LE)
idx_lower = 1:idx_LE;
x_lower = x_over_c(idx_lower);
Cp_lower = Cp(idx_lower);
ut_lower = u_t(idx_lower);

% Upper Surface (usually from LE to TE)
idx_upper = idx_LE:NPanels;
x_upper = x_over_c(idx_upper);
Cp_upper = Cp(idx_upper);
ut_upper = u_t(idx_upper);


%  PLOT TANGENTIAL VELOCITY 
figure;
plot(x_lower, ut_lower, 'r-'); % Lower Surface 
hold on;
plot(x_upper, ut_upper, 'b-'); % Upper Surface 
hold off;
xlabel('x/c');
ylabel('u_t (m/s)');
title('Tangential Component of the Velocity on the Surface');
legend('Lower Surface', 'Upper Surface', 'Location', 'best');
grid on;

%  PLOT CP DISTRIBUTION 
figure;
plot(x_lower, Cp_lower, 'r-'); % Lower Surface (Pressure Side)
hold on;
plot(x_upper, Cp_upper, 'b-'); % Upper Surface (Suction Side)
hold off;
xlabel('x/c');
ylabel('C_p');
title('C_p Distribution (Upper vs. Lower Surface)');
set(gca,'YDir','reverse'); % Flips the y-axis (standard Cp plot convention)
legend('Pressure Side (Lower)', 'Suction Side (Upper)', 'Location', 'southwest');
grid on;


%% Cp_0008_Re
A_1 = readmatrix ('cp0008_re500000.txt');
Cp_0008_500 = A_1(:,2);
X_0008_500 = A_1(:,1);

A_2 = readmatrix ('cp0008_re400000.txt');
Cp_0008_400 = A_2(:,2);
X_0008_400 = A_2(:,1);

H3 = readmatrix ('cp0008_re300000.txt');
Cp_0008_300 = H3(:,2);
X_0008_300 = H3(:,1);

figure; 

plot(X_0008_500, Cp_0008_500, 'b-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_0008_400, Cp_0008_400, 'r-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_0008_300, Cp_0008_300, 'g-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
hold off; % Releases the hold


% Set axis labels and title
xlabel('X/c (Normalized Chord Position)', 'FontSize', 12);
ylabel('C_p (Pressure Coefficient)', 'FontSize', 12);


% Add a legend to identify the lines
legend('Re = 5 x 10^5', 'Re = 4 x 10^5','Re = 3 x 10^5', 'Location', 'best');

% Invert the y-axis (Standard for Cp plots)
set(gca, 'YDir', 'reverse');

% Add a grid for easier reading
grid on;

% Value for Data Set 1 (cp0008_re500000)
Cp_0008_500_at_X0 = Cp_0008_500(1);
Cp_0008_400_at_X0 = Cp_0008_400(1);
Cp_0008_300_at_X0 = Cp_0008_300(1);


% Print the results using fprintf for clean formatting
% fprintf(' Pressure Coefficient (Cp) at Leading Edge (X/c = 0) ---\n');
% fprintf('Cp_0008 (cp0008_re500000.txt) = %.4f\n', Cp_0008_500_at_X0);
% fprintf('Cp_BD (cp_bubbledancer)    = %.4f\n', Cp_BD_at_X0);
% fprintf('------------------------------------------------------------\n');

%% Cp_BD_Re
C1 = readmatrix ('cp_bubbledancer.txt');
Cp_BD_500 = C1(:,2);
X_BD_500 = C1(:,1);

C2 = readmatrix ('cpbubbledancer_re400000.txt');
Cp_BD_400 = C2(:,2);
X_BD_400 = C2(:,1);

C3 = readmatrix ('cpbubbledancer_re300000.txt');
Cp_BD_300 = C3(:,2);
X_BD_300 = C3(:,1);

figure;

plot(X_BD_500, Cp_BD_500, 'b-', 'LineWidth', 1.5);
hold on;
plot(X_BD_400, Cp_BD_400, 'r-', 'LineWidth', 1.5);
plot(X_BD_300, Cp_BD_300, 'g-', 'LineWidth', 1.5);
hold off;


% Set axis labels and title
xlabel('X/c (Normalized Chord Position)', 'FontSize', 12);
ylabel('C_p (Pressure Coefficient)', 'FontSize', 12);


% Add a legend to identify the lines
legend('Re = 5 x 10^5', 'Re = 4 x 10^5','Re = 3 x 10^5', 'Location', 'best');

% Invert the y-axis (Standard for Cp plots)
set(gca, 'YDir', 'reverse');
grid on;

%% CP_BD_0008

B = readmatrix ('cp0008_re300000.txt');
Cp_0008 = B(:,2);
X_0008 = B(:,1);

K = readmatrix ('cpbubbledancer_re300000.txt');
Cp_BD = K(:,2);
X_BD = K(:,1);

figure;

plot(X_0008, Cp_0008, 'b-', 'LineWidth', 1.5);
hold on;
plot(X_BD, Cp_BD, 'r--', 'LineWidth', 1.5);
yline(0,'k-')
hold off;

% Set axis labels and title
xlabel('X/c (Normalized Chord Position)', 'FontSize', 12);
ylabel('C_p (Pressure Coefficient)', 'FontSize', 12);


% Add a legend to identify the lines
legend('NACA0008', 'BUBBLE DANCER', 'Location', 'best');

% Invert the y-axis (Standard for Cp plots)
set(gca, 'YDir', 'reverse');
grid on;

% 1. Find the X/c coordinate for Cp_0008 when Cp = 0
% We use linear interpolation (the 'linear' method)
% X_0008 is the independent variable, Cp_0008 is the dependent variable.
% We are interpolating the *inverse* relationship to find X where Cp is 0.
% X0_0008 = interp1(Cp_0008, X_0008, 0, 'linear');

% 2. Find the X/c coordinate for Cp_BD when Cp = 0
% X0_BD = interp1(Cp_BD, X_BD, 0, 'linear');

% 3. Print the calculated results
% fprintf('\n--- X/c Coordinate where Pressure Coefficient (Cp) = 0 ---\n');
% fprintf('X/c for Cp_0008 (cp0008_re500000) = %.4f\n', X0_0008);
% fprintf('X/c for Cp_BD (cp_bubbledancer)    = %.4f\n', X0_BD);
% fprintf('------------------------------------------------------------\n');


%% Comments
% - As it stands the Cl is correct
% - We have a 30% error on the CM but that is due to the discretization -->
%   increasing the number of panels we can fix this problem (to be done ->
%   almost done, we just have to import some files)
% - It could be interesting a plot to see how the error decays (convergence)

