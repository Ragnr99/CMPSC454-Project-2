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

% ~~~ QUESTION 3.2 ~~~
points3D = zeros(3, size(result1, 2));
for i = 1:size(result1, 2)
    % homogeneous coordinates
    p1 = [result1(:, i); 1];
    p2 = [result2(:, i); 1];
    
    % cross product matrices p1 and p2 for epipolar constraint
    crossP1 = [0, -p1(3), p1(2); p1(3), 0, -p1(1); -p1(2), p1(1), 0];
    crossP2 = [0, -p2(3), p2(2); p2(3), 0, -p2(1); -p2(2), p2(1), 0];
    
    % skew symmetric matrices applied to the camera projection matrices
    A = [crossP1 * P1; crossP2 * P2];
    
    % Singular Value Decomposition (SVD) to solve for the 3D point
    [~, ~, V] = svd(A);
    X = V(:, end); 
    X = X / X(4); % normalize
    
    points3D(:, i) = X(1:3);
end

figure;
% Plot reconstructed points
subplot(1, 2, 1);
scatter3(points3D(1, :), points3D(2, :), points3D(3, :), 'filled');
title('Reconstructed 3D Points');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; % For equal aspect ratio
grid on;

% Plot original points
subplot(1, 2, 2);
scatter3(pts3D(1, :), pts3D(2, :), pts3D(3, :), 'filled');
title('Original 3D Points');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; % For equal aspect ratio
grid on;
mse = mean(sum((points3D - pts3D).^2, 1));
sgtitle("Mean Squared Error: " + mse);