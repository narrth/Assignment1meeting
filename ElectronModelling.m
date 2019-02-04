%   ELEC 4700
%   Assignment - 1
%   Monte-Carlo Modeling of Electron Transport
%
%   Narrthanan Seevananthan

%   1 - Electron Modelling (40)


m_0 = 9.10938215e-31;   %rest mass of electrons
m_n = 0.26*m_0;         %Effective mass of electrons
k_b = 1.3806504e-23;    %Boltzmann Constant
T = 300;                %Temperature (K)
%region size is 200nm*100nm

%   A) vth (thermal velocity)

%using the following formula:
%$$vth = (k_b*T/m_n)^0.5$$
%thermal velocity was found to be 8.7053e+4 (m/s)
vth = (k_b*T/m_n)^0.5;

%   B) MFP

t_mn = 0.2e-12;     %Mean Time between collisions (ps)
mfp = vth*t_mn;     %The mean free path was found to be: 2.6449e-08 m

%using the following formula:
%$$mfp = t_mn*vth$$

%   C) Write a program that will model the random motion of the electrons
num_par = 1000;     %the number of particles in the system
time_lim = 1000;    %the max number of time steps

%randomly assign positions in the predefined space
Px = 200.*rand(1,num_par);
Py = 100.*rand(1,num_par);

%randomized direction with fixed velocity
%following pythagorean theorem we can split the velocity components
%randomly using the following formula:
% vth^2 = rand()*vx^2 + (1-rand())*vy^2
r_vel = rand (1,num_par);
v_signx = r_vel<=0.5;
v_signy = (abs(r_vel-1))<=0.5;

Vx = sqrt(r_vel.*(vth.^2));
Vy = sqrt(abs(r_vel-1).*(vth.^2));

Vx(v_signx) = -1.*Vx
Vy(v_signy) = -1.*Vy
time_step = 7.56e-15; %magic number found using vth, ydimension, and 0.01

for(t = 1:time_lim)
    Px = Px + time_step*Vx;
    Py = Py + time_step*Vx;
    
    plot(Px,Py); hold on
    xlabel('x - position');
    ylabel('y - position');
    pause(0.017)
end
%       i) 2D Plot of Particle trajectories


%       ii) Temperature Plot