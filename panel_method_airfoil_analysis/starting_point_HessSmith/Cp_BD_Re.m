C1 = readmatrix ('cp_bubbledancer.txt');
Cp_BD_500 = C1(:,2);
X_BD_500 = C1(:,1);

C2 = readmatrix ('cpbubbledancer_re400000.txt');
Cp_BD_400 = C2(:,2);
X_BD_400 = C2(:,1);

C3 = readmatrix ('cpbubbledancer_re300000.txt');
Cp_BD_300 = C3(:,2);
X_BD_300 = C3(:,1);

figure; % Opens a new figure window

plot(X_BD_500, Cp_BD_500, 'b-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_BD_400, Cp_BD_400, 'r-', 'LineWidth', 1.5);
hold on; % Keeps the first plot on the axes so the next one can be added
plot(X_BD_300, Cp_BD_300, 'g-', 'LineWidth', 1.5);
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