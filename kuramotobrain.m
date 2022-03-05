clear all
N=4;

K=@(t) sin(2*pi*rand(N,N)*t);
omega=@(t) sin(2*pi*rand(N,1)*t);
kuraN2 = @(t,theta) (kura2(t,theta,omega(t),K(t),N))';
tspan=[0,50];
cond_in=randn(N,1)*2*pi;
soluzione=ode45(kuraN2, tspan, cond_in);
delta_t=linspace(tspan(1),tspan(2),500);
theta_sol=deval(delta_t,soluzione);
% for i=1:N
%     figure
%     plot(delta_t,theta_sol(i,:));%riga i= theta_i(t)
% end

for i=1:length(delta_t)
        rer=(1/N)*sum(cos(theta_sol(:,i)));
        imr=(1/N)*sum(sin(theta_sol(:,i)));
        r(i)=sqrt(rer^2+imr^2);
end
plot(delta_t,r);
xlabel('t');
ylabel('r');
%--------------------------------------------------------------------------

function f=kura2(t,theta,omega,beta,n)
    for i=1:n
       somma=0;
       for j=1:n
           somma=somma+beta(i,j)*sin(theta(j)-theta(i));
       end 
       f(i)=omega(i)+somma/n;
    end
end
