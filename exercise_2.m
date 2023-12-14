% EXERCISE 2
% STEP 1
original_imm = imread("blister.jpg");
gray_imm = rgb2gray(original_imm);
smoothed_imm = imboxfilt(gray_imm, 15);
edge_imm = edge(smoothed_imm,'prewitt',0.007,'both');
se = strel('disk',5);
closed_imm = imclose(edge_imm,se);
stats = regionprops(closed_imm, 'BoundingBox', 'Area');
[maxArea, idx] = max([stats.Area]);
bounding_box = stats(idx).BoundingBox;
aspect_ratio = bounding_box(3) / bounding_box(4);
disp(['Aspect Ratio: ' num2str(aspect_ratio)]);
% STEP 2
gray_cropped_imm = imcrop(gray_imm, bounding_box);
adjusted_imm = imadjust(gray_cropped_imm, [0 0.68]);
hist_imm = adapthisteq(adjusted_imm);
contrast_imm = localcontrast(hist_imm, 0.2, 0.6);
adjusted_imm_2 = imadjust(contrast_imm, [0 0.6]);
hist_imm_2 = adapthisteq(adjusted_imm_2);
figure(2)
subplot(1,6,1)
imshow(gray_cropped_imm)
title("Original Cropped")
subplot(1,6,2)
imshow(adjusted_imm)
title("Adjust")
subplot(1,6,3)
imshow(hist_imm)
title("Hist Eq")
subplot(1,6,4)
imshow(contrast_imm)
title("Local Contrast")
subplot(1,6,5)
imshow(adjusted_imm_2)
title("Second Adjust")
subplot(1,6,6)
imshow(adjusted_imm_2)
title("Second Hist Eq")
% STEP 3
bin_imm = imbinarize(hist_imm_2, 0.7);
bin_imm = not(bin_imm);
se = strel("disk",3);
close_imm = imclose(bin_imm, se);
connected_elements = bwconncomp(close_imm);
max_area_index = 0;
max_area_element = 0;
for i = 1:connected_elements.NumObjects
    area_element = numel(connected_elements.PixelIdxList{i});
    if area_element > max_area_element
        max_area_element = area_element;
        max_area_label = i;
    end
end
close_background_imm = close_imm;
close_background_imm(connected_elements.PixelIdxList{max_area_label}) = 0;
small_elements_threshold = 300;
close_background_imm = bwareaopen(close_background_imm,small_elements_threshold);
min_distance_top = find(any(close_background_imm, 2), 1, 'first') - 1;
min_distance_bottom = size(close_background_imm, 1) - find(any(close_background_imm, 2), 1, 'last');
min_distance_left = find(any(close_background_imm, 1), 1, 'first') - 1;
min_distance_right = size(close_background_imm, 2) - find(any(close_background_imm, 1), 1, 'last');
if min_distance_top > min_distance_bottom
    trimmed_imm = close_background_imm(min_distance_top - min_distance_bottom + 1:end, :);
else
    trimmed_imm = close_background_imm(1:end - min_distance_bottom + min_distance_top, :);
end
if min_distance_left > min_distance_right
    trimmed_imm = trimmed_imm(:, min_distance_left - min_distance_right + 1:end);
else
    trimmed_imm = trimmed_imm(:, 1:end - min_distance_right + min_distance_left);
end
figure(3)
subplot(1,4,1)
imshow(bin_imm)
title("Binary")
subplot(1,4,2)
imshow(close_imm)
title("Closed")
subplot(1,4,3)
imshow(close_background_imm)
title("Closed no Background")
subplot(1,4,4)
imshow(trimmed_imm)
title("Trimmed")
boundaries = bwboundaries(trimmed_imm);
circularities = zeros(length(boundaries), 1);
min_length_threshold = 200;
max_length_threshold = (size(trimmed_imm,1) * 2 + size(trimmed_imm,2) * 2) * 0.5;
for k = 1:length(boundaries)
    boundary = boundaries{k};
    if length(boundary) > min_length_threshold && length(boundary) < max_length_threshold
        area = polyarea(boundary(:, 2), boundary(:, 1));
        perimeter = length(boundary);
        circularities(k) = (4 * pi * area) / (perimeter^2);
    else
        circularities(k) = 0;
    end
end
[~, max_circularity_index] = max(circularities);
max_circularity_boundary = boundaries{max_circularity_index};
min_x = min(max_circularity_boundary(:, 2));
max_x = max(max_circularity_boundary(:, 2));
min_y = min(max_circularity_boundary(:, 1));
max_y = max(max_circularity_boundary(:, 1));
pill_y = max_y - min_y; 
pill_x = max_x - min_x;
y_divisions = floor(size(trimmed_imm, 1) / (pill_y * 1.1));
x_divisions = floor(size(trimmed_imm, 2) / (pill_x * 1.1));
x_size = floor(size(trimmed_imm, 2)/x_divisions);
y_size = floor(size(trimmed_imm, 1)/y_divisions);
threshold_region_size = 400;
threshold_region_count = 3;
empty_pills = 0;
full_pills = 0;
figure(4);
imshow(trimmed_imm);
hold on;
for i = 1:x_divisions
    for j = 1:y_divisions
        x_center = floor(x_size/2) + x_size * (i - 1);
        y_center = floor(y_size/2) + y_size * (j - 1);
        
        sub_imm = trimmed_imm(y_center - floor(pill_y/2):y_center + floor(pill_y/2), x_center - floor(pill_x/2):x_center + floor(pill_x/2));
        
        sub_regions = regionprops(sub_imm, 'Area');
        num_elements = sum([sub_regions.Area] > threshold_region_size);
        
        if num_elements > threshold_region_count
            empty_pills = empty_pills + 1;
            rectangle("Position", [x_center - floor(pill_x/2), y_center - floor(pill_y/2), floor(pill_x), floor(pill_y)], 'EdgeColor', 'r')
        else
            full_pills = full_pills + 1;
            rectangle("Position", [x_center - floor(pill_x/2), y_center - floor(pill_y/2), floor(pill_x), floor(pill_y)], 'EdgeColor', 'g')
        end
        plot(x_center, y_center, "Marker", "+", "Color", "b", "MarkerSize", 5)
    end
end
hold off;
title('Subregions');
disp(['Full Pills: ' num2str(full_pills)]);
disp(['Empty Pills: ' num2str(empty_pills)]);
