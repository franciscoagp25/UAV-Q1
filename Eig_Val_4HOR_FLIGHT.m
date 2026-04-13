clear all;
clc;


g = 9.81;

% Drag matrix
Dx = 0.15;
Dy = 0.15;
Dz = 0.25;
D = diag([Dx, Dy, Dz]);

% Horizontal-flight val.s

phi_e   = 0.0;
theta_e = 0.1;     % Little value so there´s an actual roll in [rads]
psi_e   = 0.0;

vx_e = 1.0;     % in m/s
vy_e = 0.0;
vz_e = 0.0;

v_e = [vx_e; vy_e; vz_e]; %Equilibrium

%Auxiliary trigonometry values for replacement

cphi = cos(phi_e);    sphi = sin(phi_e);
cth  = cos(theta_e);  sth  = sin(theta_e);
cpsi = cos(psi_e);    spsi = sin(psi_e);

% Rotation matrix at equilibrium
R_e = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
        cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
        -sth,      sphi*cth,                   cphi*cth                  ];

% Euler-angle matrix at equilibrium
Q_e = [ 1,  sphi*tan(theta_e),  cphi*tan(theta_e);
        0,  cphi,              -sphi;
        0,  sphi/cth,           cphi/cth ];

% Matrix / Pitch dR/dphi  at equilibrium
dR_dphi = [ 0,  cphi*sth*cpsi + sphi*spsi,  -sphi*sth*cpsi + cphi*spsi;
            0,  cphi*sth*spsi - sphi*cpsi,  -sphi*sth*spsi - cphi*cpsi;
            0,  cphi*cth,                   -sphi*cth ];

% Matrix / Roll dR/dtheta  at equilibrium
dR_dtheta = [ -sth*cpsi,   sphi*cth*cpsi,   cphi*cth*cpsi;
              -sth*spsi,   sphi*cth*spsi,   cphi*cth*spsi;
              -cth,       -sphi*sth,       -cphi*sth ];

% Matrix / Yaw dR/dpsi at equilibrium
dR_dpsi = [ -cth*spsi,  -sphi*sth*spsi - cphi*cpsi,  -cphi*sth*spsi + sphi*cpsi;
             cth*cpsi,   sphi*sth*cpsi - cphi*spsi,   cphi*sth*cpsi + sphi*spsi;
             0,          0,                              0 ];

% A_p_lambda block
A_p_lambda = [dR_dphi*v_e, dR_dtheta*v_e, dR_dpsi*v_e];

% A_v_lambda block
A_v_lambda = [ 0,                    g*cos(theta_e),                  0;
              -g*cos(phi_e)*cos(theta_e),  g*sin(phi_e)*sin(theta_e), 0;
               g*sin(phi_e)*cos(theta_e),  g*cos(phi_e)*sin(theta_e), 0 ];

% S(v_e) block
S_ve = [   0,    -vz_e,   vy_e;
         vz_e,      0,   -vx_e;
        -vy_e,    vx_e,     0 ];

% Assemble A_hf
A_hf = [ zeros(3,3), A_p_lambda, R_e,         zeros(3,3);
         zeros(3,3), zeros(3,3), zeros(3,3),  Q_e;
         zeros(3,3), A_v_lambda, -D,          S_ve;
         zeros(3,3), zeros(3,3), zeros(3,3),  zeros(3,3) ];

disp('A_hf = ');
disp(A_hf);

lambda_hf = eig(A_hf);

disp('Eigenvalues of A_hf = ');
disp(lambda_hf);
