function FitToA4(fig)
    % Get current figure
    fig = gcf;
    fig.Units = 'centimeters'; 
    origSize = fig.Position(3:4); % [width, height]

    % Define A4 size limits
    a4_width = 21;  % A4 width in cm
    a4_height = 29.7; % A4 height in cm

    % Compute scaling factor to fit within A4 while maintaining aspect ratio
    scaleFactor = min(a4_width / origSize(1), a4_height / origSize(2));

    % Apply new scaled size
    newSize = origSize * scaleFactor;
    fig.Position(3:4) = newSize; % Set new width and height

    % Adjust paper size to match figure
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', newSize); % Match figure size
    set(fig, 'PaperPositionMode', 'auto'); % Prevent cropping

    % Export as high-resolution PDF
    %print(fig, 'my_figure.pdf', '-dpdf', '-r300');
return 