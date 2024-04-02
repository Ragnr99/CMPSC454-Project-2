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

% Calculated points are stored in points3D, given points are stored in pts3D

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

% ~~~ QUESTION 3.3 ~~~
% Measure 3D points on the floor
figure;
imshow(image1); 
title('Select 3 points on the floor in Image 1');
[x1, y1] = ginput(3);
close;

figure;
imshow(image2);
title('Select the same 3 points on the floor in Image 2'); 
[x2, y2] = ginput(3);
close;

% Triangulate the 3 floor points
points3D = zeros(3, 3);
for i = 1:3
    p1 = [x1(i); y1(i); 1];
    p2 = [x2(i); y2(i); 1];

    crossP1 = [0, -p1(3), p1(2); p1(3), 0, -p1(1); -p1(2), p1(1), 0];
    crossP2 = [0, -p2(3), p2(2); p2(3), 0, -p2(1); -p2(2), p2(1), 0];

    A = [crossP1 * P1; crossP2 * P2];

    [~, ~, V] = svd(A);
    X = V(:, end);
    X = X / X(4);

    points3D(:, i) = X(1:3);
end

% Fit a plane to the 3 floor points
[~, ~, V] = svd(points3D');
normalPlane = V(:, 3)'; % Last column of V is the normal vector
planePoint = mean(points3D, 2);  % Any point on the plane
planeCoefficients = [normalPlane, -dot(normalPlane, planePoint)];

% Display the floor plane equation
fprintf('plane equation: %fx + %fy + %fz + %f = 0\n', planeCoefficients);

% Define the range of x and y values
[x, y] = meshgrid(-10:0.1:10, -10:0.1:10);

% Calculate z values using the plane equation
z = (-planeCoefficients(1)*x - planeCoefficients(2)*y - planeCoefficients(4)) / planeCoefficients(3);

% Plot the surface
surf(x, y, z);

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Fitted Plane');
axis equal;
grid on;
