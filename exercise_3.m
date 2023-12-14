load loopMRI.mat
lv_contour_coordinates = cell(1, size(slice6, 4));
lv_areas = zeros(1, size(slice6, 4));
eccentricities = zeros(1, size(slice6, 4));
for i = 1:size(slice6, 4)
    % Extract the current frame
    current_frame = slice6(:, :, 1, i);
    adjusted_frame = imadjust(current_frame);
    binary_frame = adjusted_frame > 0.37 * (2^16 - 1);
    boundaries = bwboundaries(binary_frame);
    circularities = zeros(length(boundaries), 1);
    min_length_threshold = 110;
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        
        % Check if the boundary length is above the threshold
        if length(boundary) > min_length_threshold
            % Compute circularity
            area = polyarea(boundary(:, 2), boundary(:, 1));
            perimeter = length(boundary);
            circularities(k) = (4 * pi * area) / (perimeter^2);
        else
            % Assign a lower circularity value for short boundaries (e.g., 0)
            circularities(k) = 0;
        end
        [~, max_circularity_index] = max(circularities);
    end
    lv_contour_coordinates{i} = boundaries{max_circularity_index};
    lv_areas(i) = polyarea(boundaries{max_circularity_index}(:, 2) / xres, boundaries{max_circularity_index}(:, 1) / yres);
    stats = regionprops(binary_frame, 'Eccentricity');
    eccentricities(i) = stats(max_circularity_index).Eccentricity;
    
    figure;
    imshow(current_frame, []);
    hold on;
    plot(boundaries{max_circularity_index}(:, 2), boundaries{max_circularity_index}(:, 1), 'r', 'LineWidth', 2);
    title(['Frame ' num2str(i)]);
    hold off;
    %pause(0.5);
end
figure;
subplot(1,2,1)
plot(lv_areas)
title("LV Areas (mm^2)")
xlabel("Frame")
ylabel("mm^2")
subplot(1,2,2)
plot(eccentricities)
title("LV Eccentricity")
xlabel("Frame")
ylabel("Eccentricity")
save('10691904_lv_coordinates.mat', 'lv_contour_coordinates');