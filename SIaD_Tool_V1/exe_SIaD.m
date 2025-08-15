%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Q&A please mail to: luis.reyes@upc.edu

%% How to use the "Stability and Interaction assessment in
%% frequency-Domain (SIaD)" tool for modern power systems 

clear all;
clc;

%% Simulation parameters

program='Your_Project_Name.slx'; % Simulink program selection
model='Your_Project_Name'; % Main canvas name of your project (same as program usually)
Tinit=6; % Initialization time (arbitray for select the steady state)
fs=1; % Sampling frequency for FTT (specify the value in Hz for minimal freq)
delta_t=25E-6; % Fixed step time
% --- Specify the frequency vector for single-tone perturbations
fd0=unique(round(logspace(0,log10((1/delta_t)/4),3)));  % Perturbation frequencies in Hz
% --- Specify the frequency band for multi-tone or PRBS perturbations
% fd0=[0.001 6*50];  % Frequency band in Hz

%% Steady state and disturbance sources data

f0=50; % Fundamental base frequency
w=2*pi*f0; % Base angular frequency (rad/s)
Sbase=2.75E6; % Base power (VSC or System to idenfify)
Vbase=1; % Base voltage (V)
Vpeak=(Vbase/sqrt(3))*sqrt(2); % Fundamental peak voltage at the connection
Vq_ss=0; % q-component steady state voltage at bus under analysis
Vd_ss=0; % d-component steady state voltage at bus under analysis
Ipeak=Sbase/Vpeak; % Fundamental peak current through the link
Iq_ss=0; % q-component steady state current at link under analysis
Id_ss=0; % d-component steady state current at link under analysis
Vperturbation=0.03; % Percentage of nominal voltage perturbation
Iperturbation=0.03; % Percentage of nominal current perturbation

%% Optional features if a linear model is available

samples=10000; % Number of samples for linear responses comparison
jw1=1i*2*pi*(logspace(0,log10(1/delta_t),samples)); % Complex frequency vector
% Put here your linear model named as "Y_RLC1" and "Y_RLC2" or viceversa for impedance.

%% Scanning options available

% Steady-state calculation:
% 1 -> Yes, take from 1st simulation (closed loop)
% 0 -> No, I addded the d- and q-component previously (open loop)
ss_cal=1;

% Scanner settings:
% 1 -> Voltage perturbation 
% 2 -> Current perturbation 
scanner_type=1;

% 1 -> Single-tone perturbation 
% 2 -> RBS perturbation
% 3 -> Multi-tone perturbation
signal_type=3;

% 1 -> ABC scan 
% 2 -> qd0 scan 
% 3 -> pn0 scan
scanner_selector=2;

% 1 -> linear response exists
% 0 -> no linear response is available
linear=0;

% If a linear response exists, use the Optional features settings and name
% your theoretical responses as Ym_RLC or Zm_RLC for matrices of dimensions 
% of 2x2xjw1 that contain the magnitude responses of the elements and for
% the angle matrices name them as Ya_RLC or Z_RLC. SIaD automatically will
% call them inside of the master code. Put attention on the dimensions!

%% Insert here your file with your system initialization data (without run)

%------------------
% Here is your code for initizalization
Your_Initialization_Code;
%------------------

% Loading the FDs main code...
FDScanning;

% Loading the GNC Module
GNCModule;

% Loading the Modal impedance analysis Module
MAModule; % Do not use it if PRBS or Multi-tone strategy are selected!

% Loading the Phase Margin analysis Module
PMModule;

% Loading the Passivity assessment Module
PAModule;
