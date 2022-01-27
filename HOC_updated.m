%% 1-D 4nd order convection diffusion reaction equation with eps (HOC scheme)
%% Here eps or e = a
%% -alpha Uxx+bUx+cU=F(x)
%%  F(x)= -alpha *exp(((1+alpha)*(x-1))/alpha)+exp(-x)
%% exact solution
%% exp((1+alpha)*((x-1)/alpha))+e(-x) where alpha is epsilon
%% n is the number of intervals
clear all; close all; clc,
format long
alpha=1;  %Diffusion coefficient
b=1;   %Convection coefficient
x0=0;  %Lower bound
c=1;   %Reaction coefficient
xn=1;  %Upper bound
n=5;  %Total number of partition i.e total #grid=n+1
h=(xn-x0)/n;  %stepsize
A = -alpha - b^2*h^2/6/alpha + h^2*(b^2/alpha+c)/12;
B = b - h*h*b*c/12/alpha;
C =c;
x = [0:h:1]';

lowD = A/h/h - B/2/h;
d =  -2*A/h/h + c;
upD = A/h/h + B/2/h;
%%%%%%%%%#interior point=n-1 i.e unknowns%%%%%%%%%%%
I = speye(n-1,n-1);
E = sparse(2:n-1,1:n-2,1,n-1,n-1);
%M= (A*E)+(C*E')+I*B;
M= (lowD*E)+(upD*E')+I*d;
%exact solution
temp = exp((1+alpha)*(x-1)/alpha); %temp
u = temp+exp(-x);
%Source function i.s RHS
f = -alpha*temp+exp(-x);
fx = -(1+alpha)*temp - exp(-x);
fxx =  -(1+alpha)^2*temp/alpha + exp(-x);
F = f - (h*h*b/12/alpha)*fx + (h*h/12)*fxx;
F = F(2:n);
u0 = 1 + exp(-(1+alpha)/alpha);
u1 = (1+(exp(1))^-1);
F(1,1)=F(1,1)-lowD*u0;
F(n-1,1)=F(n-1,1)-upD*u1;
y=M\F;
sol = [u0;y;u1];
error = max(abs(u-sol));
%plotting the graph
fig=figure();
set(fig,'Color','white'); hold on,
plot(x,sol,'-o')
plot(x,u,'-r','LineWidth',2)
xlabel('x')
ylabel('y')
legend({'Numerical','Exact'})
grid on

% To verify convergence
% err=y-g;
% errn=norm(err,"inf")
% (log((( 5.206109099376022e-09)/(3.251867664033625e-10))))/(log(2))
