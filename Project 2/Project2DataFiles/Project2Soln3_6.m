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
% ~~~ QUESTION 3.4 ~~~
% Compute the relative rotation and translation between the two cameras
R1inv = R1'; % Inverse of a rotation matrix is its transpose
t1inv = -R1inv * t1'; % t1' converts the vector to a column vector

R2wrt1 = R2 * R1inv; % Rotation of camera 2 w.r.t. camera 1
t2wrt1 = R2 * t1inv + t2; % Translation of camera 2 w.r.t. camera 1

% Construct the essential matrix
E = R2wrt1 * skew(t2wrt1);
%disp("E: " + E)

% Compute the fundamental matrix from camera calibration parameters
Fcalib = inv(K2')' * E * inv(K1);
disp("Fcalib: " + Fcalib)


% ~~~ QUESTION 3.5 ~~~
% Measure 3D points on the floor
figure;
imshow(image1); 
title('Select 8 points on the floor in Image 1');
[x1, y1] = ginput(8);
close;

figure;
imshow(image2);
title('Select the same 8 points on the floor in Image 2'); 
[x2, y2] = ginput(8);
close;

im = image1;
im2 = image2;


figure(1); imagesc(im); axis image; drawnow;
figure(2); imagesc(im2); axis image; drawnow;


%figure(1); [x1,y1] = getpts;
figure(1); imagesc(im); axis image; hold on
for i=1:length(x1)
   h=plot(x1(i),y1(i),'*'); set(h,'Color','g','LineWidth',2);
   text(x1(i),y1(i),sprintf('%d',i));
end
hold off
drawnow;


%figure(2); imagesc(im2); axis image; drawnow;
%[x2,y2] = getpts;
figure(2); imagesc(im2); axis image; hold on
for i=1:length(x2)
   h=plot(x2(i),y2(i),'*'); set(h,'Color','g','LineWidth',2);
   text(x2(i),y2(i),sprintf('%d',i));
end
hold off
drawnow;

%do Hartley preconditioning
%savx1 = x1; savy1 = y1; savx2 = x2; savy2 = y2;
mux = mean(x1);
muy = mean(y1);
stdxy = (std(x1)+std(y1))/2;
T1 = [1 0 -mux; 0 1 -muy; 0 0 stdxy]/stdxy;
nx1 = (x1-mux)/stdxy;
ny1 = (y1-muy)/stdxy;
mux = mean(x2);
muy = mean(y2);
stdxy = (std(x2)+std(y2))/2;
T2 = [1 0 -mux; 0 1 -muy; 0 0 stdxy]/stdxy;
nx2 = (x2-mux)/stdxy;
ny2 = (y2-muy)/stdxy;

A = [];
for i=1:length(nx1);
    A(i,:) = [nx1(i)*nx2(i) nx1(i)*ny2(i) nx1(i) ny1(i)*nx2(i) ny1(i)*ny2(i) ny1(i) nx2(i) ny2(i) 1];
end
%get eigenvector associated with smallest eigenvalue of A' * A
[u,d] = eigs(A' * A,1,'SM');
F = reshape(u,3,3);

%make F rank 2
oldF = F;
[U,D,V] = svd(F);
D(3,3) = 0;
F = U * D * V';

%unnormalize F to undo the effects of Hartley preconditioning
F = T2' * F * T1;

colors =  'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';
%overlay epipolar lines on im2
L = F * [x1' ; y1'; ones(size(x1'))];
[nr,nc,nb] = size(im2);
figure(2); clf; imagesc(im2); axis image;
hold on; plot(x2,y2,'*'); hold off
for i=1:length(L)
    a = L(1,i); b = L(2,i); c=L(3,i);
    if (abs(a) > (abs(b)))
       ylo=0; yhi=nr; 
       xlo = (-b * ylo - c) / a;
       xhi = (-b * yhi - c) / a;
       hold on
       h=plot([xlo; xhi],[ylo; yhi]);
       set(h,'Color',colors(i),'LineWidth',2);
       hold off
       drawnow;
    else
       xlo=0; xhi=nc; 
       ylo = (-a * xlo - c) / b;
       yhi = (-a * xhi - c) / b;
       hold on
       h=plot([xlo; xhi],[ylo; yhi],'b');
       set(h,'Color',colors(i),'LineWidth',2);
       hold off
       drawnow;
    end
end


%overlay epipolar lines on im1
L = ([x2' ; y2'; ones(size(x2'))]' * F)' ;
[nr,nc,nb] = size(im);
figure(1); clf; imagesc(im); axis image;
hold on; plot(x1,y1,'*'); hold off
for i=1:length(L)
    a = L(1,i); b = L(2,i); c=L(3,i);
    if (abs(a) > (abs(b)))
       ylo=0; yhi=nr; 
       xlo = (-b * ylo - c) / a;
       xhi = (-b * yhi - c) / a;
       hold on
       h=plot([xlo; xhi],[ylo; yhi],'b');
       set(h,'Color',colors(i),'LineWidth',2);
       hold off
       drawnow;
    else
       xlo=0; xhi=nc; 
       ylo = (-a * xlo - c) / b;
       yhi = (-a * xhi - c) / b;
       hold on
       h=plot([xlo; xhi],[ylo; yhi],'b');
       set(h,'Color',colors(i),'LineWidth',2);
       hold off
       drawnow;
    end
end

disp("F.2: " + F)


% ~~~ QUESTION 3.6.1 ~~~
% Initialize variables to accumulate squared distances
sedCalib = 0;
sedReal = 0;

% SED for FCalib
for i = 1:numel(x1)
    % epipolar lines
    lineIn2From1 = Fcalib * [x1(i); y1(i); 1];
    lineIn1From2 = Fcalib' * [x2(i); y2(i); 1];
    
    % squared distances
    squaredDist2From1 = (lineIn2From1(1) * x2(i) + lineIn2From1(2) * y2(i) + lineIn2From1(3))^2 / (lineIn2From1(1)^2 + lineIn2From1(2)^2);
    squaredDist1From2 = (lineIn1From2(1) * x1(i) + lineIn1From2(2) * y1(i) + lineIn1From2(3))^2 / (lineIn1From2(1)^2 + lineIn1From2(2)^2);
    
    sedCalib = sedCalib + squaredDist2From1 + squaredDist1From2;
end

% mean SED error for Fcalib
sedCalib = sedCalib / (2 * numel(x1));

% SED for F
for i = 1:numel(x1)
    % epipolar lines
    lineIn2From1 = F * [x1(i); y1(i); 1];
    lineIn1From2 = F' * [x2(i); y2(i); 1];
    
    % compute squared distances
    squaredDist2From1 = (lineIn2From1(1) * x2(i) + lineIn2From1(2) * y2(i) + lineIn2From1(3))^2 / (lineIn2From1(1)^2 + lineIn2From1(2)^2);
    squaredDist1From2 = (lineIn1From2(1) * x1(i) + lineIn1From2(2) * y1(i) + lineIn1From2(3))^2 / (lineIn1From2(1)^2 + lineIn1From2(2)^2);
    
    sedReal = sedReal + squaredDist2From1 + squaredDist1From2;
end

% mean SED error for F
sedReal = sedReal / (2 * numel(x1));

fprintf('SED error for Fcalib: %f\n', sedCalib);
fprintf('SED error for F: %f\n', sedReal);

function S = skew(v)
    S = [0    -v(3)  v(2);
         v(3)   0   -v(1);
        -v(2)  v(1)    0];
end