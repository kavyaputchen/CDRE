%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-D 2th order convection diffusion reaction equation with eps exponential scheme (EXP-CDRE)  
%% Here eps or e = a                                                                           
%% -aUxx+bUx+cU=F(x)                                                                           
%% exact solution                                                                               
%% exp((1+a)*((x-1)/a))+e(-x) where a is epsilon                                                
%% F(x)=  -a*((exp(((1+a)*((x)-1))/a))+exp((-x)))                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%format long
a=1e-3;
b=1;
x0=0; %left boundary
c=1;
xn=1;  %right boundary
n=40;  %number of interval
h=(xn-x0)/n;  %step size

x = [x0:h:xn]';

P=((b*h)/(2*a));
alpha=a*P*coth(P);
A=-(alpha/(h*h))-(b/(2*h));
B=(2*alpha/(h*h))+c;
C=-(alpha/(h*h))+(b/(2*h));
I = speye(n-1,n-1);

F = -a*((exp(((1+a)*((x)-1))/a))+exp((-x)));

%Forming A matrix
E = sparse(2:n-1,1:n-2,ones(n-2,1),n-1,n-1);

M= (A*E)+(C*E')+I*B;

u0 = (1+exp(-(1+a)/a)); % influence from the left boundary
un = (1+(exp(1))^-1); % influence from the right boundary

%%%%%%%%%%%%%Exact solution%%%%%%%%%%%%%
g = ((exp(((1+a)*((x)-1))/a))+exp((-x)));

f = F(2:n);
f(1,1)=f(1,1)-(A*u0);
f(n-1,1)=f(n-1,1)-(C*un);

u=M\f;
sol =[u0;u;un]; 

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
%(log(((8.944390452680917e-06)/(4.453866530973460e-06))))/(log(2))
toc
