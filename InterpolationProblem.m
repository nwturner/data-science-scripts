function [C] = InterpolationProblem(x,y,Num)
%C=InterpolationProblem(x,y,Num) finds the coefficients of the interpolation function
%based on a set of data.
% x and y are vectors representing a set of data. Num is the number of
% data points based on which we wish to find our interpolation function.
N = Num;
x0=linspace(min(x),max(x),Num); % New x vector based on Num data points
y0=zeros(1,length(x0)); % Initializing new y vector based on Num data points
step=(length(x)-1)/(Num-1);
for i=1:length(x0)
y0(i)=y((i-1)*step+1);
end
A = zeros(N,N); % Initializing matrix A corresponding to the interpolation function
for i=1:N % i is the row index
for j=1:N % j is the column index
A(i,j)=x0(i)^(N-j); % Assigning entries to the matrix A
end
end
C=A\y0'; % Solving the equation A*C = y0' for the interpolation function coefficients
end