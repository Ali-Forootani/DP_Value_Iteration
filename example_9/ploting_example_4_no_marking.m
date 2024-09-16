%% Plotting for Example 1

% Set the figure size to be larger
figure('Position', [100, 100, 1000, 800]);

% First subplot
subplot(4,1,1)
hold on
plot(frequency_decision_BD(:,1)/st_size(1), 'LineWidth', 2)
%plot(frequency_decision_net(:,1)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
ylabel('C_1', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
%legend('Booking Case', 'No Booking Case', 'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Second subplot
subplot(4,1,2)
hold on
plot(frequency_decision_BD(:,2)/st_size(1), 'LineWidth', 2)
%plot(frequency_decision_net(:,2)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
ylabel('C_2', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Third subplot
subplot(4,1,3)
hold on
plot(frequency_decision_BD(:,3)/st_size(1), 'LineWidth', 2)
%plot(frequency_decision_net(:,3)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('C_3', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Optional Fourth subplot
subplot(4,1,4)
hold on
plot(frequency_decision_BD(:,4)/st_size(1), 'LineWidth', 2)
% plot(frequency_decision_net(:,4)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('C_4', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% General settings
set(gcf, 'Color', 'w'); % Set the background to white

% Apply font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

% Add title to the figure
sgtitle('Frequency distribution of each action', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');


%%

%% Data Extraction for Revenue Plot

% Extract the data you want to plot
data1 = squeeze(Revenue(10, 5, :));
data2 = squeeze(Revenue(30, 5, :));
data3 = squeeze(Revenue(50, 5, :));
data4 = squeeze(Revenue(70, 5, :));

% Extract title values for each subplot
title_values1 = Revenue(10, 1:4);
title_values2 = Revenue(30, 1:4);
title_values3 = Revenue(50, 1:4);
title_values4 = Revenue(70, 1:4);

% Set the figure size and background
figure('Position', [100, 100, 1200, 800]);
set(gcf, 'Color', 'w');

% First subplot
subplot(2, 2, 1);
plot(data1, 'LineWidth', 2);
xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('(%.0f, %.0f, %.0f, %.0f)', title_values1), 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
grid on;

% Second subplot
subplot(2, 2, 2);
plot(data2, 'LineWidth', 2);
xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('(%.0f, %.0f, %.0f, %.0f)', title_values2), 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
grid on;

% Third subplot
subplot(2, 2, 3);
plot(data3, 'LineWidth', 2);
xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('(%.0f, %.0f, %.0f, %.0f)', title_values3), 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
grid on;

% Fourth subplot
subplot(2, 2, 4);
plot(data4, 'LineWidth', 2);
xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('(%.0f, %.0f, %.0f, %.0f)', title_values4), 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
grid on;

% Apply consistent font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');


%% Data Extraction Loop Version

% Extract the data and title values for each subplot
data = cell(1, 4);
titles = cell(1, 4);
for i = 1:4
    data{i} = squeeze(Revenue(10*i, 5, :));
    titles{i} = Revenue(10*i, 1:4);
end

% Set up the figure
figure('Position', [100, 100, 1200, 800]);
set(gcf, 'Color', 'w');

% Plot and annotate each subplot
for i = 1:4
    subplot(2, 2, i);
    plotAndAnnotate(1:length(data{i}), data{i}, titles{i}, i);
end

% Apply consistent font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

% Define the function to plot and annotate
function plotAndAnnotate(x, y, title_values, subplotIndex)
    plot(y, 'LineWidth', 2);
    
    % Label and title settings
    xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold'); % Smaller font size for y-axis label
    title(sprintf('(%.0f, %.0f, %.0f, %.0f)', title_values), 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    
    grid on;
    
    % Highlight the first value with a different marker
    hold on;
    firstValue = y(1);
    plot(1, firstValue, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Red circle for first value with larger size
    
    % Calculate offset to place the annotation slightly above and to the right of the point
    verticalOffset = -0.09 * range(y); % Vertical offset
    horizontalOffset = 5.2; % Horizontal offset in terms of x-axis units
    
    % Annotate to the right of the first value
    text(1 + horizontalOffset, firstValue + verticalOffset, sprintf(' %.2f', firstValue), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 18, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', 'k'); % Black text for visibility
    
    hold off;

    % Adjust Y-axis limits if necessary
    ylim([min(y) - 0.1*range(y), max(y) + 0.2*range(y)]); % Increased limit for visibility
end
