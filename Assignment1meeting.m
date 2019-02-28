%   ELEC 4700
%   Assignment - 1
%   Monte-Carlo Modeling of Electron Transport
%
%   Narrthanan Seevananthan

%   1 - Electron Modelling (40)
clear

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
num_par = 10;     %the number of particles in the system
time_lim = 100;    %the max number of time steps

%randomly assign positions in the predefined space
Px = 200e-9.*rand(1,num_par);
Py = 100e-9.*rand(1,num_par);

%randomized direction with fixed velocity
%following pythagorean theorem we can split the velocity components
%randomly using the following formula:
% vth^2 = rand()*vx^2 + (1-rand())*vy^2
r_vel = rand (1,num_par);
v_signx = rand (1,num_par)<=0.5;
v_signy = rand (1,num_par)<=0.5;

Vx = sqrt(r_vel.*(vth.^2));
Vy = sqrt(abs(r_vel-1).*(vth.^2));

Vx(v_signx) = -1.*Vx(v_signx);
Vy(v_signy) = -1.*Vy(v_signy);

%       i) 2D Plot of Particle trajectories
time_step = 0.1*(100e-9)/vth; %magic number found using vth, ydimension, and 0.01
cmap = hsv(num_par);    %color map to plot different color lines


P_setx = zeros(num_par,time_lim);
P_sety = zeros(num_par,time_lim);

sc_temp = (1:time_lim)*0;

for t = 1:time_lim
    P_setx(:,t) = Px.';
    P_sety(:,t) = Py.';
    

    sc_temp(t) = sum(Vx.^2 + (Vy).^2)*m_n/(k_b*num_par);

    
    for i = 1:num_par
     
        if (0<Px(i)) && (Px(i)<(200e-9)) && (t>1)
            figure(1);
            plot(P_setx(i,(t-1):t),P_sety(i,(t-1):t),'Color',cmap(i,:));     
            hold on
            
        end
    end
    
    Px = Px + time_step*Vx;
    Py = Py + time_step*Vy;
    
    Px(Px>200e-9) = 0;
    Px(Px<0) = 200e-9;
    
    Vy(Py>100e-9) = -1.*Vy(Py>100e-9);
    Vy(Py<0) = -1.*Vy(Py<0);
    
    title('2-D Electron Trajectories')
    xlabel('x - position');
    ylabel('y - position');
    xlim([0 200e-9]);       %x dimensions in plot
    ylim([0 100e-9]);       %y dimensions in plot
    
    %       ii) Temperature Plot
    figure(2);
    plot(1:time_lim,sc_temp);
    title('semiconductor temperature')
    xlabel('Time');
    ylabel('Temperature');
    
    pause(0.0017)
end



