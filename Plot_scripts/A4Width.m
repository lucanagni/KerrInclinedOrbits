function A4Width(fig)
    fig.Units = 'centimeters'; 
    origSize = fig.Position(3); % [width]

    % Define A4 size limits
    a4_width = 20.9;  % A4 width in cm

    % Compute scaling factor to fit within A4 while maintaining aspect ratio
    scaleFactor = a4_width / origSize(1);

    % Apply new scaled size
    newSize = origSize * scaleFactor;
    fig.Position(3) = newSize; % Set new width
    %fig.Position(4) = fig.Position(4)*scaleFactor;

    % Adjust paper size to match figure
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', fig.Position(3:4)); % Match figure size
    set(fig, 'PaperPositionMode', 'auto'); % Prevent cropping

    % Export as high-resolution PDF
    %print(fig, 'my_figure.pdf', '-dpdf', '-r300');
return
