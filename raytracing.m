close all
clear all
clc

%% Rotations
Rx = @(a) [
    1,      0,       0;
    0, cos(a), -sin(a);
    0, sin(a),  cos(a);
];
Ry = @(a) [
    cos(a), 0,  sin(a);
    0,      1,       0;
   -sin(a), 0,  cos(a);
];
Rz = @(a) [
    cos(a), -sin(a)  0;
    sin(a),  cos(a)  0;
    0,       0,      1;
];

Rtot = @(a,b,c) Rz(c)*Ry(b)*Rx(a);

%% PPBIO ray tracing

% Surfaces
% bottom
Rb = eye(3);
Ob = [0;0;0];

% side1
Rs1 = Rx(34/180*pi);
Os1 = [0;-11.4;1.5];


% top
Rt = eye(3);
Ot = [0;0;3];

% side2
Rs2 = Rx(-34/180*pi);
Os2 = [0;11.4;1.5];

% lens
Rl = eye(3);
Ol = [0;-11.4;-3];
fx = 12.3;
fy = 12.3;

% rays
N = 100;
[k0,P0] = gen_rays(eye(3),Ol - [0;0;5],1.8,1.8,N,0);

% trace
[k1,P1] = lens(k0,P0,Rl,Ol,fx,fy);
[k2,P2] = refract(k1,P1,Rb,Ob,1,1.45);
[k3,P3] = reflect(k2,P2,Rs1,Os1);
[k4,P4] = reflect(k3,P3,Rb,Ob);
[k5,P5] = reflect(k4,P4,Rt,Ot);

rays = [];
for i =1:N
    rays = [rays, P0(:,i), P1(:,i), P2(:,i), P3(:,i), P4(:,i), P5(:,i),NaN(3,1)];
end

plot3(rays(1,:),rays(2,:),rays(3,:)); hold on; axis equal

Pi = 0.5*[[1;1;0],[-1;1;0],[-1;-1;0],[1;-1;0],[1;1;0]];

cpb = Rb*(diag([20,26.72,1])*Pi) + Ob;
cps1 = Rs1*(diag([20,5.359,1])*Pi) + Os1;
cpt = Rt*(diag([20,17.84,1])*Pi) + Ot;
cps2 = Rs2*(diag([20,5.359,1])*Pi) + Os2;
cpl = Rl*(diag([10,10,1])*Pi) + Ol;

planes = [cpb,NaN(3,1),cps1,NaN(3,1),cpt,NaN(3,1),cps2,NaN(3,1),cpl];

plot3(planes(1,:),planes(2,:),planes(3,:))
