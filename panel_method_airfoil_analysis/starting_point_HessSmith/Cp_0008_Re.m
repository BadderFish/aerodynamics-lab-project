A_1 = readmatrix ('cp0008_re500000.txt');
Cp_0008_500 = A_1(:,2);
X_0008_500 = A_1(:,1);

A_2 = readmatrix ('cp0008_re400000.txt');
Cp_0008_400 = A_2(:,2);
X_0008_400 = A_2(:,1);

H3 = readmatrix ('cp0008_re300000.txt');
Cp_0008_300 = H3(:,2);
X_0008_300 = H3(:,1);

figure; % Opens a new figure window

plot(X_0008_500, Cp_0008_500, 'b-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_0008_400, Cp_0008_400, 'r-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_0008_300, Cp_0008_300, 'g-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
hold off; % Releases the hold

% --- Customizing the Plot ---

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
fprintf('--- Pressure Coefficient (Cp) at Leading Edge (X/c = 0) ---\n');
fprintf('Cp_0008 (cp0008_re500000) = %.4f\n', Cp_0008_at_X0);
fprintf('Cp_BD (cp_bubbledancer)    = %.4f\n', Cp_BD_at_X0);
fprintf('------------------------------------------------------------\n');