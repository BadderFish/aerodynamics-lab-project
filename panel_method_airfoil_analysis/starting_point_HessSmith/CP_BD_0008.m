B = readmatrix ('cp0008_re500000.txt');
Cp_0008 = B(:,2);
X_0008 = B(:,1);

K = readmatrix ('cp_bubbledancer.txt');
Cp_BD = K(:,2);
X_BD = K(:,1);
% --- Plotting the Data ---

figure; % Opens a new figure window

plot(X_0008, Cp_0008, 'b-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_BD, Cp_BD, 'r--', 'LineWidth', 1.5);
yline(0,'k-')
hold off; % Releases the hold

% --- Customizing the Plot ---

% Set axis labels and title
xlabel('X/c (Normalized Chord Position)', 'FontSize', 12);
ylabel('C_p (Pressure Coefficient)', 'FontSize', 12);


% Add a legend to identify the lines
legend('NACA0008', 'BUBBLE DANCER', 'Location', 'best');

% Invert the y-axis (Standard for Cp plots)
set(gca, 'YDir', 'reverse');

% Add a grid for easier reading
grid on;

% 1. Find the X/c coordinate for Cp_0008 when Cp = 0
% We use linear interpolation (the 'linear' method)
% X_0008 is the independent variable, Cp_0008 is the dependent variable.
% We are interpolating the *inverse* relationship to find X where Cp is 0.
X0_0008 = interp1(Cp_0008, X_0008, 0, 'linear');

% 2. Find the X/c coordinate for Cp_BD when Cp = 0
X0_BD = interp1(Cp_BD, X_BD, 0, 'linear');

% 3. Print the calculated results
fprintf('\n--- X/c Coordinate where Pressure Coefficient (Cp) = 0 ---\n');
fprintf('X/c for Cp_0008 (cp0008_re500000) = %.4f\n', X0_0008);
fprintf('X/c for Cp_BD (cp_bubbledancer)    = %.4f\n', X0_BD);
fprintf('------------------------------------------------------------\n');