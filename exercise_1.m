% EXERCISE 1
% STEP 1
data = load("hand.mat");
original_imm = data.J;
%equalized_imm = adapthisteq(original_hand);
%smoothed_imm = medfilt2(original_imm);
%smoothed_imm = imboxfilt(original_imm, 3);
A = 1.7;
kernel_size =11;
blurred = imgaussfilt(original_imm, kernel_size);
high_boost = A * double(original_imm) - double(blurred);
high_boost = uint8(max(0, min(255, high_boost)));
adjusted_imm = imadjust(high_boost,[0 0.55]);
smoothed_imm = imboxfilt(adjusted_imm, 3);
%

%
edge_imm = smoothed_imm;
% Display the transformed image
figure(1)
subplot(2,3,1)
imshow(original_imm)
title("Original")
subplot(2,3,2)
imshow(high_boost)
title("High Boost of smoothed")
subplot(2,3,3)
imshow(adjusted_imm)
title("Adjusted")
subplot(2,3,4)
imshow(smoothed_imm)
title("Smoothed")
subplot(2,3,5)
%
subplot(2,3,6)
imshow(edge_imm)
title("Final Edge Enhancement")
% STEP 2 - Robinson and Roberts
compass_masks = [-1 -2 -1; 0 0 0; 1 2 1]; 
edge_robinson = imfilter(adjusted_imm, compass_masks, 'conv');
threshold = 0.045; 
edge_roberts = edge(edge_imm, 'Roberts',threshold,'nothinning');
figure(2)
subplot(1,2,1)
imshow(edge_robinson)
title("Robinson")
subplot(1,2,2)
imshow(edge_roberts)
title("Roberts")
% STEP 3 - Morphological Operators
% Robinson
% Robert
se_1 = strel("disk",3);
se_2 = strel("line",3,0);
robert_close = bwareaopen(edge_roberts,15);
robert_close = bwfill(robert_close,"holes",4);
robert_open = imopen(robert_close, se_2);
figure(4)
subplot(1,3,1)
imshow(edge_roberts)
title("Original Roberts")
subplot(1,3,2)
imshow(robert_close)
title("Close Robert")
subplot(1,3,3)
imshow(robert_open)
title("Open Robert")
% STEP 4 - Visualize Results
robert_boundaries = bwboundaries(robert_open);
figure(6)
imshow(original_imm);
hold on;
for k = 1:length(robert_boundaries)
   boundary = robert_boundaries{k};
   plot(boundary(:,2), boundary(:,1), 'r')
end
hold off
