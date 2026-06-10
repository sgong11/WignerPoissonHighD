% --- Configuration ---
% Assuming f_ht is a cell array containing data for all time steps
num_snapshots = length(f_ht); 
Nx = 512;
% If you have an actual array of time values (e.g., t_array), uncomment and use it.
% Otherwise, the code will just display the snapshot index.
% t_array = linspace(0, 10, num_snapshots); % Example time array

% Generate the meshgrid once outside the loop
[X, V] = meshgrid(x{1}, v{1});

% --- Video Setup ---
video_filename = '0513.mp4';
vObj = VideoWriter(video_filename, 'MPEG-4'); % Use 'MPEG-4' for .mp4 or 'Motion JPEG AVI' for .avi
vObj.FrameRate = 15; % Adjust playback speed (frames per second)
open(vObj);

% Create and format the figure before the loop
fig = figure;

% --- Loop over time snapshots ---
for k = 1:num_snapshots
    
    % 1. Extract data for the current snapshot (using k instead of end)
    % A = ht.evaluate_slice(f_ht{k}, 1, 256, 1:Nx, 1:Nx);
    A4 = zeros(Nx, Nx); 
    for j = 1:Nx
        A4(:, j) = ht.evaluate_fiber(f_ht{k}, 1, j, 1:Nx,j);
    end
    
    % 2. Plot the contour
    contourf(X, V, A4', 20, 'LineColor', 'none');
    colorbar; % Optional: helps track amplitude changes over time
    
    % 3. Update the title with the current time
    % If using a time array, change to: title(sprintf('Time = %.3f', t_array(k)));
    title(sprintf('Time Snapshot: %f', t(k))); 
    
    xlabel('X');
    ylabel('V');
    
    % Force MATLAB to draw the updated plot
    drawnow;
    
    % 4. Capture the figure as a frame and write it to the video file
    frame = getframe(fig);
    writeVideo(vObj, frame);
    
end

% --- Clean up ---
close(vObj);
disp(['Movie successfully saved as ', video_filename]);