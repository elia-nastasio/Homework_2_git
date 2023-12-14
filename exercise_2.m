% EXERCISE 2
% STEP 1
original_imm = imread("blister.jpg");
gray_imm = rgb2gray(original_imm);
smoothed_imm = imboxfilt(gray_imm, 15);
edge_imm = edge(smoothed_imm,'prewitt',0.007,'both');
se = strel('disk',5);
close_imm = imclose(edge_imm,se);
stats = regionprops(close_imm, 'BoundingBox', 'Area');
[maxArea, idx] = max([stats.Area]);
bounding_box = stats(idx).BoundingBox;
aspect_ratio = bounding_box(3) / bounding_box(4);
disp(['Aspect Ratio: ' num2str(aspect_ratio)]);
%{
figure(1)
subplot(1,4,1)
imshow(gray_imm)
title("Original")
subplot(1,4,2)
imshow(edge_imm)
title("Edges")
subplot(1,4,3)
imshow(close_imm)
title("Closed")
subplot(1,4,4)
imshow(gray_imm)
hold on
rectangle("position",bounding_box,"EdgeColor",'r')
hold off
title("Boudning box")
%}
% STEP 2
gray_cropped_imm = imcrop(gray_imm, bounding_box);
adjusted_imm = imadjust(gray_cropped_imm, [0 0.68]);
hist_imm = adapthisteq(adjusted_imm);
contrast_imm = localcontrast(hist_imm, 0.2, 0.6);
adjusted_imm_2 = imadjust(contrast_imm, [0 0.6]);
hist_imm_2 = adapthisteq(adjusted_imm_2);

%{
sharpened_image = imsharpen(adjusted_image);
contrast_imm = localcontrast(adjusted_image, 0.4, 0.7);
%}
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
%{
subplot(1,6,3)
imshow(adjusted_image)
title("Adjusted")
subplot(1,6,4)
imshow(adapt_equalized_imm)
title("Adapt Hist")
subplot(1,6,5)
imshow(sharpened_image)
title("Sharpened")
subplot(1,6,6)
imshow(contrast_imm)
title("Local Contrast Processing")
%}
% STEP 3
%T = graythresh(adjusted_imm);
bin_imm = imbinarize(hist_imm_2, 0.7);
bin_imm = not(bin_imm);
se = strel("disk",3);
close_imm = imclose(bin_imm, se);
close_back_imm = close_imm;
background_region = regionprops(close_imm, 'Area', 'PixelIdxList');
[~, idx_bckground] = max([background_region.Area]);
biggest_region = zeros(size(close_start_imm));
biggest_region(background_region(idx).PixelIdxList) = 1;
biggest_region_inverse = ~biggest_region;
close_imm(biggest_region_inverse) = 0;
figure(3)
subplot(1,3,1)
imshow(bin_imm)
title("Binary")
subplot(1,3,2)
imshow(close_start_imm)
title("Closed")
subplot(1,3,3)
imshow(close_imm)
title("Closed with no Background")
%
boundaries = bwboundaries(close_imm);
circularities = zeros(length(boundaries), 1);
min_length_threshold = 200;
max_length_threshold = (size(close_imm,1) * 2 + size(close_imm,2) * 2) * 0.9;
% Loop through boundaries
for k = 1:length(boundaries)
    boundary = boundaries{k};
    if length(boundary) > min_length_threshold & length(boundary) < max_length_threshold
        % Compute circularity
        area = polyarea(boundary(:, 2), boundary(:, 1));
        perimeter = length(boundary);
        circularities(k) = (4 * pi * area) / (perimeter^2);
    else
        circularities(k) = 0;
    end
end

% Find the indices of the first two elements with max circularity
[sorted_circularities, sorted_indices] = sort(circularities, 'descend');
N = 1;
max_circularity_indices = sorted_indices(1:N);

% Display the original image with the identified objects
figure;
imshow(close_imm);
hold on;
%bounding_box_pills = zeros(200, 1);
filled_pills = 0;
empty_pills = 0;
for i = 1:N
    % Plot the boundary of the object
    max_circularity_boundary = boundaries{max_circularity_indices(i)};
    min_x = min(max_circularity_boundary(:, 2));
    max_x = max(max_circularity_boundary(:, 2));
    min_y = min(max_circularity_boundary(:, 1));
    max_y = max(max_circularity_boundary(:, 1));
    pill_y = max_y - min_y; 
    pill_x = max_x - min_x;
    y_divisions = floor(size(close_imm, 1) / (pill_y*1.1));
    x_divisions = floor(size(close_imm, 2) / (pill_x*1.1));
    rect_width = floor(size(close_imm, 2) / x_divisions);
    rect_height = floor(size(close_imm, 1) / y_divisions);

    for j = 1:x_divisions
        for k = 1:y_divisions
            % Calculate the position of the rectangle with respect to (0, 0) of the image
            rect_x_start = (j - 1) * rect_width + 1;
            rect_y_start = (k - 1) * rect_height + 1;
            
            % Ensure that the indices stay within bounds
            rect_x_end = rect_x_start + rect_width;
            rect_y_end = rect_y_start + rect_height;
            
            % Extract the sub-image within the rectangle
            sub_image = close_imm(rect_y_start:rect_y_end, rect_x_start:rect_x_end);

            % Count the number of elements in the sub-image
            stats_3 = regionprops(sub_image, 'Area');

            % Set a threshold for the number of elements
            threshold = 3000; % Adjust this threshold as needed
            num_elements = sum([stats_3.Area] > threshold);
            % Color the rectangle based on the number of elements
            if num_elements > 4
                empty_pills = empty_pills + 1;
                rectangle('Position', [rect_x_start, rect_y_start, rect_width, rect_height], 'EdgeColor', 'r');
            else
                filled_pills = filled_pills + 1;
                rectangle('Position', [rect_x_start, rect_y_start, rect_width, rect_height], 'EdgeColor', 'g');
            end
            % Plot the rectangle
            %rectangle('Position', [rect_x, rect_y, rect_width, rect_height], 'EdgeColor', 'b');
        end
    end

    plot(max_circularity_boundary(:, 2), max_circularity_boundary(:, 1), 'r', 'LineWidth', 2);
    %rectangle('Position', [min_x, min_y, max_x - min_x, max_y - min_y], 'EdgeColor', 'b');
end

hold off;
title('Objects with Top 2 Circularities and Satisfying Length Threshold');
disp(['Filled Pills: ' num2str(filled_pills)]);
disp(['Empty Pills: ' num2str(empty_pills)]);
%{
cc = bwconncomp(thin_imm);
stats = regionprops(cc, 'Area', 'Eccentricity', 'BoundingBox');
unusedPillAreaThreshold = 500;
usedPillEccentricityThreshold = 0.9;
unusedPillMask = false(size(thin_imm));
usedPillMask = false(size(thin_imm));
for i = 1:length(stats)
    if stats(i).Area > unusedPillAreaThreshold && stats(i).Eccentricity < usedPillEccentricityThreshold
        unusedPillMask(cc.PixelIdxList{i}) = true;
    elseif stats(i).Eccentricity >= usedPillEccentricityThreshold
        usedPillMask(cc.PixelIdxList{i}) = true;
    end
end
figure(4);
subplot(1, 3, 1);
imshow(thin_imm);
title('Thin');

subplot(1, 3, 2);
imshow(unusedPillMask);
title('Unused Pill Mask');

subplot(1, 3, 3);
imshow(usedPillMask);
title('Used Pill Mask');
%}