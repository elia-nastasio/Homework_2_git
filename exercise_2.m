% EXERCISE 2
% STEP 1
original_imm = imread("blister.jpg");
gray_imm = rgb2gray(original_imm);
smoothed_imm = imboxfilt(gray_imm, 25);
edge_imm = edge(smoothed_imm,'prewitt',0.007,'both');
%figure(1)
%imshow(edge_imm)
%title("Edge image");
se = strel('disk',5);
close_imm = imclose(edge_imm,se);
%figure(2)
%imshow(close_imm)
%title('Closed image')
stats = regionprops(close_imm, 'BoundingBox', 'Area');
[maxArea, idx] = max([stats.Area]);
bounding_box = stats(idx).BoundingBox;
aspect_ratio = bounding_box(3) / bounding_box(4);
disp(['Aspect Ratio: ' num2str(aspect_ratio)]);
% STEP 2
cropped_imm = imcrop(original_imm, bounding_box);
gray_cropped_imm = rgb2gray(cropped_imm);
figure(3)
subplot(1,2,1)
imshow(gray_cropped_imm)
title("Original Image")
equalized_imm = adapthisteq(gray_cropped_imm);
subplot(1,2,2)
imshow(equalized_imm)
title("adapthisteq")
% STEP 3
T = graythresh(equalized_imm);
A = imbinarize(equalized_imm, T);
A = not(A);
rowPercentages = mean(A, 2);
colPercentages = mean(A, 1);
threshold = 0.7;
rowsToZero = rowPercentages > threshold;
colsToZero = colPercentages > threshold;
A(rowsToZero, :) = 0;
A(:, colsToZero) = 0;
se_2 = strel("disk",11);
se_3 = strel("disk",11);
cropped_open = imopen(A, se_2);
cropped_close = imclose(cropped_open, se_3);
thin_imm = bwmorph(cropped_close,"thin");
%dilated_shrink = imdilate(thin_imm,se_2);
figure(5)
subplot(2,2,1)
imshow(A)
title("Original binary")
subplot(2,2,2)
imshow(cropped_open)
title("Open")
subplot(2,2,3)
imshow(cropped_close)
title("Closed")
subplot(2,2,4)
imshow(thin_imm)
title("Thin");
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
figure(6);
subplot(1, 3, 1);
imshow(thin_imm);
title('Thin');

subplot(1, 3, 2);
imshow(unusedPillMask);
title('Unused Pill Mask');

subplot(1, 3, 3);
imshow(usedPillMask);
title('Used Pill Mask');
