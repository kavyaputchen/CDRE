%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  1-D 2th order convection diffusion reaction equation with eps (CD-CDRE)
%% Here eps or epsilon is a
%% -aUxx+bUx+cU=F(x)
%% F(x)= -a*exp(((1+a)*(x-1))/a)+exp(-x)
%% Exact solution
%% exp((1+a)*((x-1)/a))+e(-x) where a is epsilon
%% n is the number of intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
format long
a=1;
b=1;
x0=0;   %% Left boundary
c=1;
xn=1;   %% Right boundary
n=40;    %% No of intervals
h=(xn-x0)/n; %% step size
x=[x0:h:xn]';

A=(-a/(h*h))+(-b/(2*h));
B=((2*a)/(h*h))+c;
C=(-a/(h*h))+(b/(2*h));

I = speye(n-1,n-1);

F=-a*((exp(((1+a)*((x)-1))/a))+exp((-x)));

%Forming A matrix
E = sparse(2:n-1,1:n-2,ones(n-2,1),n-1,n-1);
M= (A*E)+(C*E')+I*B;

%%%%%%%%%%%%%%%%%%%%%exact solution%%%%%%%%%%%%%%%%%%%%%%%%

g=((exp(((1+a)*((x)-1))/a))+exp((-x)));

%Forming B matrix

U0 = (1+exp(-(1+a)/a)); %%left boundary condition
Un = (1+(exp(1))^-1); %%Right boundary condition

f=F(2:n);
f(1,1)=F(1,1)-(A*U0);
f(n-1,1)=(F(n-1,1))-(C*Un);

u=M\f;
sol =[U0;u;Un]; 

%plotting the graph
fig=figure(); hold on
set(fig,'Color','white')
plot(x,sol,'-o','LineWidth',1,'MarkerSize',8)
plot(x,g,'LineWidth',2)
xlabel('x')
ylabel('y')
legend({'Exponential','Exact'})
grid on

% To verify convergence
err=sol-g;
errn=norm(err,"inf")
%(log(((1.574523575724740e-04)/(2.296746763841284e-05))))/(log(2))
toc