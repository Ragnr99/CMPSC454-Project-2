clear;
load("mocapPoints3D.mat");

% convert 3D world coordinates to homogeneous coordinates
worldHomog = [pts3D; ones(1, size(pts3D, 2))];

% parameters_V1
load("Parameters_V1.mat");
image1 = imread('im1corrected.jpg');
K1 = Parameters.Kmat; 
R1 = Parameters.Rmat; 
t1 = Parameters.position;

% parameters_V2
load("Parameters_V2.mat");
image2 = imread('im2corrected.jpg');
K2 = Parameters.Kmat; 
R2 = Parameters.Rmat; 
t2 = Parameters.position;

% ~~~ QUESTION 3.1 ~~~

% convert 3D world coordinates to homogeneous coordinates
worldHomog = [pts3D; ones(1, size(pts3D, 2))];

% P = K * [R|t]
P1 = K1 * [R1, -R1 * t1'];

% project from 3D world coordinates to 2D image coordinates
imageHomog1 = P1 * worldHomog;

% Normalize homogeneous coordinates
x1 = imageHomog1(1, :) ./ imageHomog1(3, :);
y1 = imageHomog1(2, :) ./ imageHomog1(3, :);

result1 = [x1; y1];

% P = K * [R|t]
P2 = K2 * [R2, -R2 * t2'];

% project from 3D world coordinates to 2D image coordinates
imageHomog2 = P2 * worldHomog;

% normalize homogeneous coordinates
x2 = imageHomog2(1, :) ./ imageHomog2(3, :);
y2 = imageHomog2(2, :) ./ imageHomog2(3, :);

result2 = [x2; y2];

figure;
% first image
subplot(1, 2, 1);
imshow(image1);
hold on;
plot(result1(1, :), result1(2, :), 'ro', 'MarkerSize', 2);
hold off;
title('Camera 1');

% second image
subplot(1, 2, 2);
imshow(image2);
hold on;
plot(result2(1, :), result2(2, :), 'ro', 'MarkerSize', 2);
hold off;
title('Camera 2');