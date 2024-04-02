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
% ~~~ QUESTION 3.4 ~~~
% Compute the relative rotation and translation between the two cameras
R1inv = R1'; % Inverse of a rotation matrix is its transpose
t1inv = -R1inv * t1'; % t1' converts the vector to a column vector

R2wrt1 = R2 * R1inv; % Rotation of camera 2 w.r.t. camera 1
t2wrt1 = R2 * t1inv + t2; % Translation of camera 2 w.r.t. camera 1

% Construct the essential matrix
E = R2wrt1 * skew(t2wrt1);
disp(E)

% Compute the fundamental matrix from camera calibration parameters
Fcalib = inv(K2')' * E * inv(K1);
disp(Fcalib)

function S = skew(v)
    S = [0    -v(3)  v(2);
         v(3)   0   -v(1);
        -v(2)  v(1)    0];
end