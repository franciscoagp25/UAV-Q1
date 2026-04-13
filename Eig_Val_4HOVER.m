clear all;
clc;


g = 9.81;

% Drag matrix
Dx = 0.15;
Dy = 0.15;
Dz = 0.25;
D = diag([Dx, Dy, Dz]);

% Hover  A matrix
A = [ zeros(3,3), zeros(3,3), eye(3),    zeros(3,3);
      zeros(3,3), zeros(3,3), zeros(3,3), eye(3);
      zeros(3,3), [0 g 0; -g 0 0; 0 0 0], -D, zeros(3,3);
      zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3) ];

disp('A = ');
disp(A);

lambda = eig(A);

disp('Eigenvalues of A = ');
disp(lambda);
