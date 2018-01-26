%fdtd1d_step1.m
% Basic FDTD Engine
%initialize matlab/octave
close all;
clc;
clear all;

%units MKS
meters = 1;
centimeters = 1e-2*meters;
millimeters = 1e-3*meters;
inches = 2.54 * centimeters;
feet = 12 * inches;

seconds = 1;

hertz = 1 / seconds;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;

%constants
c0 = 29979458 * meters/seconds; %speed of light
e0 = 8.8541878176e-12 * 1 / meters; %permitivity
u0 = 1.2566370614e-6 * 1 / meters; %permeability

%Open a figure
figure('color', 'w')

%%%%%%%% Dashboard %%%%%%%%%%%%%%
%FDTD parameters
dz = 0.006 * meters; % grid resolution in the z direction
Nz = 200; % 200 points in the grid
dt = 1e-11 * seconds; %time step
STEPS = 1000; %1000 time steps

%%%%%% Buid device on grid %%%%%
%initialize material to free space
ER = ones(1, Nz);%relative permeability air
UR = ones(1, Nz);%relative permitivity

%%%% Initialize FDTD parameters %%%%%%
%Cimpute update coefficients
mEy = (c0*dt)./ER; 
mHx = (c0*dt)./UR; 

%Initialize fields
Ey = zeros(1, Nz);
Hx = zeros(1, Nz);


%%%% Perform FDTD Analysis %%%%%
%%% Main Loop
for T = 1 : STEPS
  % Update H from E
  for nz = 1: Nz-1
      Hz(nz) = Hx(nz) + mHx(nz)*( Ey(nz+1) - Ey(nz))/dz;
  endfor
  Hz(Nz) = Hx(Nz) + mHx(Nz)*( 0 - Ey(Nz))/dz;
  %Update E from H
  Ey(1) = Ey(1) + mEy(1)*(Hx(1) - 0)/dz;
  for nz = 2 : Nz
    Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/dz;
  endfor
  
  %Show status
  if ~mod(T,10)
    %show fields
    %draw1d(ER, Ey, Hz, dz);
    hor_axis = dz:dz:Nz*dz;
    plot(hor_axis, Ey, hor_axis, Hx);
    xlim([dz Nz*dz]);
    xlabel('z');
    title(['Field at step ' num2str(T) ' of ' num2str(STEPS)]);
    drawnow
  endif
  
endfor


