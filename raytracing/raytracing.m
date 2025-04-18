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
Os1 = [0;-11.4;1.5]*1e-3;


% top
Rt = eye(3);
Ot = [0;0;3]*1e-3;

% side2
Rs2 = Rx(-34/180*pi);
Os2 = [0;11.4;1.5]*1e-3;

% collimation lens
Rl1 = eye(3);
Ol1 = [0;-11.302;-10]*1e-3;
fx1 = 5.1*1e-3;
fy1 = 5.1*1e-3;

% foccusing lens
Rl2 = eye(3);
Ol2 = [0;-11.302;-3]*1e-3;
fx2 = Inf;
fy2 = 12.49*1e-3;

% rays
% [k0,P0] = gen_rays(eye(3),Ol - [0;0;5],1.8,1.8,100,0);
tx0 = 20*pi/180;
ty0 = 20*pi/180;
[k0,P0] = gen_gaussian_beam(eye(3), Ol1 - [0;0;5e-3],tx0,ty0,0.8e-6,4);

% trace
[k1,P1] = lens(k0,P0,Rl1,Ol1,fx1,fy1);
[k2,P2] = lens(k1,P1,Rl2,Ol2,fx2,fy2);
[k3,P3] = refract(k2,P2,Rb,Ob,1,1.45);
[k4,P4] = reflect(k3,P3,Rs1,Os1);
[k5,P5] = reflect(k4,P4,Rb,Ob);
[k6,P6] = reflect(k5,P5,Rt,Ot);
[k7,P7] = reflect(k6,P6,Rb,Ob);
[k8,P8] = reflect(k7,P7,Rs2,Os2);
[k9,P9] = reflect(k8,P8,Rb, Ob);


figure
plot(P6(1,:),P6(2,:),'r.');
% axis equal

%% Plot PPBIO
rays = [];
for i = 1:1:length(k0)
    rays = [rays, P0(:,i), P1(:,i), P2(:,i), P3(:,i), P4(:,i), P5(:,i), P6(:,i), P7(:,i), P8(:,i), P9(:,i) ,NaN(3,1)];
end


figure
plot3(rays(1,:),rays(2,:),rays(3,:),'r'); hold on; axis equal

Pi = 0.5*[[1;1;0],[-1;1;0],[-1;-1;0],[1;-1;0],[1;1;0]];

cpb = Rb*(1e-3*diag([20,26.72,1])*Pi) + Ob;
cps1 = Rs1*(1e-3*diag([20,5.359,1])*Pi) + Os1;
cpt = Rt*(1e-3*diag([20,17.84,1])*Pi) + Ot;
cps2 = Rs2*(1e-3*diag([20,5.359,1])*Pi) + Os2;
cpl1 = Rl1*(1e-3*diag([10,10,1])*Pi) + Ol1;
cpl2 = Rl2*(1e-3*diag([10,10,1])*Pi) + Ol2;

planes = [cpb,NaN(3,1),cps1,NaN(3,1),cpt,NaN(3,1),cps2,NaN(3,1),cpl1,NaN(3,1),cpl2];

plot3(planes(1,:),planes(2,:),planes(3,:))

figure
plot(rays(2,:),rays(3,:),'r'); hold on;
plot(planes(2,:),planes(3,:)); axis equal

