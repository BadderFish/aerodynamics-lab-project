%% test Hess Smith

clc
close all
clear 

addpath mat_functions

%% Input

U_inf = 1;      % Far-field velocity [m/s]
AoA = -3;       % Angle of attack (you have to take it negative because it is the free stream angle of attack)
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
[xx, ind_sort] = sort(x_over_c);

figure;
plot(xx, u_t(ind_sort), 'o');
xlabel('x/c');
ylabel('u_t (m/s)');
title('Tangential component of the velocity on the surface');
grid on;

figure;
plot(xx, Cp(ind_sort), 'o');
xlabel('x/c');
ylabel('C_p');
title('C_p distribution');
set(gca,'YDir','reverse'); % We print -Cp
grid on;


%% Comments
% - As it stands the Cl is correct
% - We have a 30% error on the CM but that is due to the discretization -->
%   increasing the number of panels we can fix this problem (to be done ->
%   almost done, we just have to import some files)
% - It could be interesting a plot to see how the error decays (convergence)


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


% --- PLOT TANGENTIAL VELOCITY ---
figure;
plot(x_lower, ut_lower, 'b-'); % Lower Surface in blue
hold on;
plot(x_upper, ut_upper, 'r-'); % Upper Surface in red
hold off;
xlabel('x/c');
ylabel('u_t (m/s)');
title('Tangential Component of the Velocity on the Surface');
legend('Lower Surface', 'Upper Surface', 'Location', 'best');
grid on;

% --- PLOT CP DISTRIBUTION ---
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


