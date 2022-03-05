clear all
N=50;
K=2;

omega=randn(N,1);
iniz=rand(N,1)*2*pi;
%iniz=zeros(N,1);

kuraN = @(t,theta) (kura(t,theta,omega,K,N))';
tspan=[0,20];

soluzione=ode45(kuraN, tspan, iniz);
%[T, theta_sol] = ode45(kuraN, tspan , omega);
delta_t=linspace(tspan(1),tspan(2),200);
theta_sol=deval(delta_t,soluzione);

 for i=1:N
      plot(delta_t,theta_sol(i,:),"b");%riga i= theta_i(t)
      hold on;
      axis([19,20,-7,7])
      xlabel('t') 
      ylabel('\theta_i(t)')
 end
 %faccio il grafico delle traiettorie delle fasi degli oscillatori nel 
 %tempo
 
 
 
 
 %-------------------------------------------------------------------------
 function f=kura(t,theta,omega,beta,n)
   for i=1:n
       somma=0;
       for j=1:n
           somma=somma+sin(theta(j)-theta(i));
       end 
       f(i)=omega(i)+(beta/n)*somma;
   end
 end 