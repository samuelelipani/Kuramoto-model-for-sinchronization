%in questo programma desidero rappresentare nello spazio delle fasi
%l'evoluzione della fase media del parametro d'ordine pertanto dovrò
%innanzitutto trovare le soluzioni del sistema di equazioni differenziali
%di kuramoto, scrivere un vettore delle derivate di psi in funzione di un
%certo intervallino tstep e graficarlo in funzione di psi.
clear all
N=5;
K=2;
%utilizzo un coupling che mi garantisca una sincronizzazione della
%popolazione di oscillatori
omega=randn(N,1);
iniz=rand(N,1)*2*pi;
%le frequenze distribuite gaussianamente e le condizioni inziali
%equidistribuite sulla circonferenza unitaria
kuraN = @(t,theta) (kura(t,theta,omega,K,N))';
tspan=[0,100];
soluzione=ode45(kuraN, tspan, iniz);
delta_t=linspace(tspan(1),tspan(2),100000);
theta_sol=deval(delta_t,soluzione);

psi=[];
dpsidt=[];
tstep=0.0001;
for k=1:length(delta_t)
    r=(1/N*sum(exp(i*theta_sol(:,k))));
    psi(k) = 2*atan(imag(r)/(real(r))); 
end
dpsidt=diff(psi)/tstep;
dpsidt(1,100000)=0;
plot(psi,dpsidt);
xlabel('psi') 
ylabel('dpsi/dt')
legend('traiettoria del sistema nello spazio delle fasi')





%--------------------------------------------------------------------------
function f=kura(t,theta,omega,beta,n)
   for i=1:n
       somma=0;
       for j=1:n
           somma=somma+sin(theta(j)-theta(i));
       end 
       f(i)=omega(i)+(beta/n)*somma;
   end
end