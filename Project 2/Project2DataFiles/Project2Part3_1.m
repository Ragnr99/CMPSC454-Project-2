clear;
load("mocapPoints3D.mat");
%load("Parameters_V1.mat");
load("Parameters_V2.mat");
%image = imread('im1corrected.jpg');
image = imread('im2corrected.jpg');

% convert 3D world coordinates to homogeneous coordinates
worldHomog = [pts3D; ones(1, size(pts3D, 2))];

K = Parameters.Kmat; 
R = Parameters.Rmat; 
t = Parameters.position;

% P = K * [R|t]
Pmat = K * [R, -R * t'];

% project from 3D world coordinates to 2D image coordinates
imageHomog = Pmat * worldHomog;

% normalize homogeneous coordinates
x = imageHomog(1, :) ./ imageHomog(3, :);
y = imageHomog(2, :) ./ imageHomog(3, :);

result = [x; y];

imshow(image);
hold on;
plot(result(1, :), result(2, :), 'ro', 'MarkerSize', 5); % Plot red points with a marker size of 5
hold off;