%fdtd1d_step2.m
% Add Gaussian source to FDTD Engine
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
%Source Parameters
fmax = 5.0 * gigahertz;

%Grid Parameters
nmax = 1; %max refraction index
NLAM = 10; % grid resolution of ten points
NBUFZ = [100 100]; %buffer zone, 100 points before and after.

%%%%% Compute opimized grid %%%%%%
%Nominal resolution
lam0 = c0/fmax; %wavelength at max freq;
dz = lam0 / nmax / NLAM;

%Compute grid size
Nz = sum(NBUFZ) + 3; %number of point on the grid

% Comupute Grid Axis
za = [0:Nz-1] * dz; %grid axis

%%%%%% Buid device on grid %%%%%
%initialize material to free space
ER = ones(1, Nz);%relative permeability air
UR = ones(1, Nz);%relative permitivity

%%% Compute the source
%compute time step (dt)
nbc = sqrt(UR(1)*ER(1));
dt = nbc*dz/(2*c0); % one new cell in two time clicks

%compute source parameters
tau = 0.5/fmax;
t0 = 5*tau;

%Compute the number of time steps
tprop =nmax*Nz*dz/c0; %time it takes to traverse the grid
t = 2*t0 + 3*tprop;
STEPS = ceil(t/dt); %rounded total number of steps

%Compute the source
t = [0:STEPS-1]*dt; %time axis (reused variable t)
nz_src = round(Nz/2); %position of the source
Esrc = exp(-((t- t0)/tau).^2);

%%%% Initialize FDTD parameters %%%%%%
%Cimpute update coefficients
mEy = (c0*dt)./ER; 
mHx = (c0*dt)./UR; 

%Initialize fields
Ey = zeros(1, Nz);
Hx = zeros(1, Nz);


%%%% Perform FDTD Analysis %%%%%
%%% Main Loop
for T = 1 :  STEPS
  % Update H from E
  for nz = 1: Nz-1
    Hx(nz) = Hx(nz) + mHx(nz)*( Ey(nz+1) - Ey(nz))/dz;
  end
  Hx(Nz) = Hx(Nz) + mHx(Nz)*( 0 - Ey(Nz))/dz;
  %Update E from H
  Ey(1) = Ey(1) + mEy(1)*(Hx(1) - 0)/dz;
  for nz = 2 : Nz
    Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/dz;
  end
  
  %inject source
  Ey(nz_src) = Ey(nz_src) + Esrc(T);
  
  %Show status
  if ~mod(T,5)
    %show fields
    %draw1d(ER, Ey, Hz, dz);
    hor_axis = dz:dz:Nz*dz;
    plot(za, Ey, 'b', za, Hx, 'r', linewidth=2);
    xlim([dz Nz*dz]);
    xlabel('z');
    title(['Field at step ' num2str(T) ' of ' num2str(STEPS)]);
    drawnow
  endif
endfor