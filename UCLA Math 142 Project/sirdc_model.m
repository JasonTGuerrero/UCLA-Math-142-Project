%Code for simulating an SIRDC model
clear all;
close all;

totTime = 12; %total simulation time
dt = 0.03; %time step
tvec = 0:dt:totTime;
tsteps = length(tvec); %total number of time steps

%initialize vectors for susceptibles, infected, and recovered
S = zeros(tsteps,1); 
I = zeros(tsteps,1);
R = zeros(tsteps,1);
D = zeros(tsteps,1);
C = zeros(tsteps,1);

shi = 1; %re-infection probability
theta = 0.1; %incubation average
delta = 0.04; %death rate
beta = 0.24; %contact rate (you choose this!)
gamma = 0.96; %recovery rate (you choose this!)
N = 10000; %total population

%Initialize susceptibles and infected. We choose I(1) equal to some significant nonzero value for visibility.
S(1) = 9500; 
I(1) = 500;
R(1) = 0;
D(1) = 0;
C(1) = 0;

for i=1:tsteps-1
    S(i+1) = S(i) + dt*((-beta*S(i)*I(i)/N) + shi*C(i));
    I(i+1) = I(i) + dt*(beta*S(i)*I(i)/N - gamma*I(i));
    R(i+1) = R(i) + dt*(gamma*I(i) - theta*R(i));
    D(i+1) = D(i) + dt*(delta*theta*R(i));
    C(i+1) = C(i) + dt*((1-delta)*theta*R(i) - shi*C(i));
end

%Plot S, I, R, D and C over time. Feel free to edit this figure to your liking.
figure(1)
plot(tvec,S,'r',tvec,I,'b',tvec,R,'g',tvec, D, 'k',tvec, C, 'c','linewidth',2)
xlabel('time')
ylabel('population')
