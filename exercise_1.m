% EXERCISE 1
% STEP 1
data = load("hand.mat");
original_imm = data.J;
%equalized_imm = adapthisteq(original_hand);
smoothed_imm = medfilt2(original_imm);
%adjusted_imm = imadjust(smoothed_imm,[0.5 0.55], [0.7 0.75]);
adjusted_imm = adapthisteq(smoothed_imm);
% Display the transformed image
figure(1)
enhanced_imm = adjusted_imm;
imshowpair(original_imm,adjusted_imm,"montage")
title('Original and Edge Enhancement');
[Gmag,Gdir] = imgradient(adjusted_imm,"prewitt");
figure(3)
imshowpair(Gmag,Gdir,"montage");
% STEP 2
%robinson_adjusted = adapthisteq(edge_robinson);
%edge_imm = imadjust(smoothed_imm);
compass_masks = [-1 -2 -1; 0 0 0; 1 2 1]; 
edge_robinson = imfilter(adjusted_imm, compass_masks, 'conv');
threshold = 0.006; 
sigma = 10;
edge_canny = edge(adjusted_imm,'canny',threshold, sigma);
figure(2)
subplot(1,2,1)
imshow(edge_robinson)
title("Robinson")
subplot(1,2,2)
imshow(edge_canny)
title("Canny")
%subplot(1,4,4)
%imshow(robinson_adjusted)
%title("Robinson hist adjusted")
