% Uses Euler method to solve chemical reaction formula given as a diffential equation and plots the concentration of each
% chemical component as a function of time for various reaction rate constants. The equation solved here is d[H20]/dt = k[OH][H], 
% where k is the reaction rate constant.


a=0; % start time
b=0.2; % end time
dx=0.01; 
kf=[20,30,40]; % array of kf values
H2Oin = 0.2; % initial concentration of H2O
[t1,H2O_1,OH_1,H_1] = ODESolverth(@Eulerth,@RHSth,a,b,H2Oin,dx,kf(2));
[t2,H2O_2,OH_2,H_2] = ODESolverth(@Eulerth,@RHSth,a,b,H2Oin,dx,kf(1));
[t3,H2O_3,OH_3,H_3] = ODESolverth(@Eulerth,@RHSth,a,b,H2Oin,dx,kf(3));

% Printing final [H2O] for different kf values
fprintf('Final [H2O] for kf=%d: %f\n',kf(2),H2O_1(end,1))
fprintf('Final [H2O] for kf=%d: %f\n',kf(1),H2O_2(end,1))
fprintf('Final [H2O] for kf=%d: %f\n',kf(3),H2O_3(end,1))

% Plotting [H2O]
figure(1)
plot(t1,H2O_1,t2,H2O_2,t3,H2O_3)
xlabel('Time (s)'); ylabel('H20 concentration');
title('[H20] vs Time'); legend('kf=30','kf=20','kf=40')
% Plotting [OH]
figure(2)
plot(t2,OH_1,t2,OH_2,t3,OH_3)
xlabel('Time (s)'); ylabel('OH concentration');
title('[OH] vs Time'); legend('kf=30','kf=20','kf=40')
% Plotting [H]
figure(3)
plot(t3,H_1,t2,H_2,t3,H_3)
xlabel('Time (s)'); ylabel('H concentration');
title('[H] vs Time'); legend('kf=30','kf=20','kf=40')

function [X,Y] = Eulerth(X0,Y0,DX,RHS,OH,H,kf)
%Euler method
F = RHS(OH,H,kf);
Y = Y0 + DX * F;
X = X0 + DX;
end

function [X,Y,arrayOH,arrayH] = ODESolverth(Integrator,RHS,A,B,YA,DX,kf)
%Numerical ODE solver
OH = 0.5;
H = 0.4;
NI = int64(((B-A)/DX)+1); % Number of integration steps
X=zeros(NI,1);
Y=zeros(NI,1);
% Initial condition
X(1)=A;
Y(1)=YA;
arrayOH(1) = OH;
arrayH(1) = H;
for i=1:NI-1 % Now we can perform NI integration steps
    [X(i+1),Y(i+1)] = Integrator(X(i),Y(i),DX,RHS,OH,H,kf);
    OH = OH - ( Y(i+1)-Y(i) );
    arrayOH(i+1) = OH;
    H = H - ( Y(i+1)-Y(i) );
    arrayH(i+1) = H;
end
end

function [F] = RHSth(OH,H,kf)
%RHS of d[H20]/dt = k[OH][H]
%   Right hand side
F = kf*OH*H;
end
