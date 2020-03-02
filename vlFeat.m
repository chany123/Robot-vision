close all;
clear all;

% Load images.
im0 = imread('im0.png');
im1 = imread('im1.png');
imS0 = single(rgb2gray(im0));
imS1 = single(rgb2gray(im1));
% Calculate SIFT descriptors. 
[f0, d0] = vl_sift(imS0);
[f1, d1] = vl_sift(imS1);
% Calculate descriptor correspondences. 
[matches01, scores01] = vl_ubcmatch(d0, d1);
% Sort distance of two points by descending. 
[drop, perm] = sort(scores01, 'descend');
matches01 = matches01(:, perm);
% scores  = scores(perm);
% Drawing the descriptor correspondences of im0 and im1.
draw(im0, im1, f0, f1, matches01, imS0, 'im0 and im1');

% Load images.
im1L = imread('im1L.png');
imS1L = single(rgb2gray(im1L));
% Calculate SIFT descriptors. 
[f1L, d1L] = vl_sift(imS1L);
% Calculate descriptor correspondences. 
[matches01L, scores01L] = vl_ubcmatch(d0, d1L);
% Sort distance of two points by descending. 
[drop, perm] = sort(scores01L, 'descend');
matches01L = matches01L(:, perm);
% Drawing the descriptor correspondences of im0 and im1L.
draw(im0, im1L, f0, f1L, matches01L, imS0, 'im0 and im1L');

% Calculate and show the translation and overlaid images of im0 and im1.
figure; subplot(1, 2, 1);
bestV01 = transAndOver(f0, f1, matches01, im0, im1, 'Output of part 3.2. 01.');
% Calculate and show the translation and overlaid images of im0 and im1L.
subplot(1, 2, 2);
bestV01L = transAndOver(f0, f1L, matches01L, im0, im1L, 'Output of part 3.2. 01L.');

% Load im1E_crop1 and im1E_crop2image.
im1E_crop1 = imread('im1E_crop1.png');
im1E_crop2 = imread('im1E_crop2.png');
% im1E_crop1 = padarray(im1E_crop1, [1599 2631], 'post'); 
imSE_crop1 = single(rgb2gray(im1E_crop1));
imSE_crop2 = single(rgb2gray(im1E_crop2));

% Calculate SIFT descriptors. 
[fc1, dc1] = vl_sift(imSE_crop1);
[fc2, dc2] = vl_sift(imSE_crop2);
% Calculate descriptor correspondences. 
[matches0c1, scores0c1] = vl_ubcmatch(d0, dc1, 5);
[matches0c2, scores0c2] = vl_ubcmatch(d0, dc2, 5);

% Calculate and show the translation, rotation and overlaid images of im0 and im1E_crop1.
figure; subplot(1, 2, 1);
[bestV0c1, bestR0c1] = transRotAndOver(f0, fc1, matches0c1, imS0, imSE_crop1, 'Output of part 3.3. 01.');
subplot(1, 2, 2);
[bestV0c2, bestR0c2] = transRotAndOver(f0, fc2, matches0c2, imS0, imSE_crop2, 'Output of part 3.3. 02.');

% The function of drawing the descriptor correspondences of two images.
function draw(imageA, imageB, fa, fb, matches, ima, name) 
    figure; 
    subplot(2, 1, 1);
    imagesc(cat(2, imageA, imageB)); title(['Original image of ', name, '.']);
    axis image off ;
    subplot(2, 1, 2); imagesc(cat(2, imageA, imageB));
    hold on
    xa = fa(1, matches(1, :));
    xb = fb(1, matches(2, :)) + size(ima, 2);
    ya = fa(2, matches(1, :));
    yb = fb(2, matches(2, :));
    h = line([xa ; xb], [ya ; yb]);
    set(h, 'color', 'b', 'linewidth', 1);
    vl_plotframe(fa(:,matches(1, :)));
    fb(1, :) = fb(1, :) + size(ima, 2);
    vl_plotframe(fb(:, matches(2, :)));
    title(['Descriptor correspondences of ', name, '.']);
    axis image off ;
end

% The function of calculate and show the translation and overlaid images of two images.
function bestV = transAndOver(fout, fin, matches, imout, imin, name) 
    % Get the match data points.
    dataA = fout(1:2, matches(1, :));
    dataB = fin(1:2, matches(2, :));
    number = size(dataA, 2);
    % Set the repeat times.
    iter = 500;
    % Set the allowed error.
    sigma = 10;
    pretotal = -1;
    for i = 1:iter
        % Get a random point.
        idx = randperm(number, 1); 
        sampleA = dataA(:, idx);
        sampleB = dataB(:, idx); 

        % Calculate the translation vector from the input image to the output image.
        transV = [(sampleA(1, :) - sampleB(1,:)); (sampleA(2, :) - sampleB(2,:))];
        transB = dataB + transV;

        % Calculate the distance between the transformed image and the original image.
        dis = sqrt((dataA(1, :) - transB(1,:)).^2 + (dataA(2, :) - transB(2,:)).^2);
        % Total the number of points whose error is lower than the set allowed error.
        total = sum(dis < sigma);              

        if total > pretotal            
            pretotal = total;
            % Find the best translation vector.
            bestV = transV;          
        end  
    end
    % Set the transformation matrice. 
    T = [1 0 0; 0 1 0; bestV(1) bestV(2) 1];
    transT = affine2d(T);
    % Translation and Overlaid images. 
    outputView = imref2d(size(imout));
    transR = imwarp(imin(:,:,1), transT, 'OutputView', outputView);
    transG = imwarp(imin(:,:,2), transT, 'OutputView', outputView);
    transB = imwarp(imin(:,:,3), transT, 'OutputView', outputView);
    transI(:,:,1) = transR;
    transI(:,:,2) = transG;
    transI(:,:,3) = transB;
    % Show the picture with the different colour channel.
    transI = imfuse(imout, transI, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 0 2]);
    imshow(transI); title(name);
    axis image off ;
end

% The function of calculate and show the translation, rotation and overlaid images of two images.
function [bestV, bestR] = transRotAndOver(fout, fin, matches, imout, imin, name) 
    % Get the match data points.
    dataA = fout(:, matches(1, :));
    dataB = fin(:, matches(2, :));
    number = size(dataA, 2);
    % Set the repeat times.
    iter = 45;
    % Set the allowed error.
    sigma = 30;
    pretotal = -1;

    for i = 1:iter
        % Get a random point.
        idx = randperm(number, 1); 
        sampleA = dataA(1:2, idx);
        sampleB = dataB(1:2, idx); 

        % Calculate the rotation angle. 
        transB = dataB(1:2, :);
        transR = dataA(4, idx) - dataB(4, idx);
        transM = [cos(transR) -sin(transR); sin(transR) cos(transR)];
        % The rotated image. 
        for j = 1:size(transB, 2)
            transB(:, j) = transM * transB(:, j);
        end
    
        % Calculate the translation vector from the input image to the output image.
        transV = [(sampleA(1) - sampleB(1)); (sampleA(2) - sampleB(2))];
        transV = transV + sampleB - transM * sampleB;
        transB = transB + transV;

        % Calculate the distance between the transformed image and the original image.
        dis = sqrt((dataA(1, :) - transB(1,:)).^2 + (dataA(2, :) - transB(2,:)).^2);
        % Total the number of points whose error is lower than the set allowed error.
        total = sum(dis < sigma);              

        if total > pretotal            
            pretotal = total;
            % Find the best translation vector.
            bestV = transV;    
            % Find the best rotation radian.
            bestR = transR;
        end  
    end

    C = createOverlaidImage(imout, imin, bestV, bestR);
    imshow(C); title(name);
    axis image off ;
end

