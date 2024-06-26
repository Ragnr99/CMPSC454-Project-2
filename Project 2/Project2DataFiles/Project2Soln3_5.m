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


for j=1:3
    for i=1:3
        fprintf('%10g ',10000*F(j,i));
    end
    fprintf('\n');
end



