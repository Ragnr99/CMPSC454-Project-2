% Load camera parameters
load('Parameters_V1.mat');

% Load 3D mocap points
load('mocapPoints3D.mat');

K = data_v1.Parameters.Kmat; 
R = data_v1.Parameters.Rmat; 
t = data_v1.Parameters.position;
P = K * [R, t'];

pinhole = zeros(3, num_points);
point = zeros(4,1)

% Perform 3D to 2D projection for all points
for i = 1:num_points
    for(j = 1:3) 
        point(j,1) = (pts3D(j, i));
        point(4,1) = 1;
    end
    result = P * point;  
    pinhole(:, i) = result
end

% Load the image
image = imread('im1corrected.jpg');
imshow(image)



