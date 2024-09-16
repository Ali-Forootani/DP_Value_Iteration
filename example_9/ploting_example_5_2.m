%% Data Extraction
% Extract the data and title values for each subplot
data = cell(1, 5);
titles = cell(1, 5);
for i = 1:5
    data{i} = squeeze(Revenue(i, 5, :));
    titles{i} = Revenue(i, 1:4);
end
% Add the 5th dataset and title values
data{5} = squeeze(Revenue(5, 5, :));
titles{5} = Revenue(5, 1:4);

% Set up the figure
figure('Position', [100, 100, 1200, 1500]); % Adjust height for 5 subplots
set(gcf, 'Color', 'w');

% Plot and annotate each subplot in a 5x1 layout
for i = 1:5
    subplot(5, 1, i); % Change to 5x1 subplot layout
    plotAndAnnotate(1:length(data{i}), data{i}, titles{i}, i);
end

% Apply consistent font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

% Define the function to plot and annotate
function plotAndAnnotate(x, y, title_values, subplotIndex)
    plot(y, 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 2);
    
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
    verticalOffset = 0.05 * range(y); % Vertical offset
    horizontalOffset = 5; % Horizontal offset in terms of x-axis units
    
    % Annotate to the right of the first value
    text(1 + horizontalOffset, firstValue + verticalOffset, sprintf(' %.2f', firstValue), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 18, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', 'k'); % Black text for visibility
    
    hold off;

    % Adjust Y-axis limits if necessary
    ylim([min(y) - 0.1*range(y), max(y) + 0.2*range(y)]); % Increased limit for visibility
end
