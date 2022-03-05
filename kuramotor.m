clear 
clf
%Vogliamo rappresentare l'evoluzione temporale in funzione della costante
%di accoppiamento del parametro d'ordine. Calcolato numericamente e
%confrontandolo con quello analitico che viene fuori da una distribuzione
%delle frequenze degli oscillatori lorentziana 
N=1000; 
K=[1,2];
%mi serve avere due valori della costante di accoppiamento in modo poi da
%plottare due curve per i due diversi parametri d'ordine 
omega=randn(N,1);
iniz=[zeros(N,1),rand(N,1)*2*pi];
%iniz=[rand(N,1)*2*pi,rand(N,1)*2*pi];
%provo poi a mettere al posto di zeros i numeri random e vedo che si
%attesta a zero sempre perchè il parametro d'ordine parte già quasi nullo
%essendo il nucleo di oscillatori non sincronizzati incoerenti così in
%partenza!
%nel calcolo del parametro d'ordine quando k è minore di kc posso anche
%fare partire tutti gli oscillatori da fasi iniziali nulle poichè non
%essendo possibile la sincronizzazione inizieranno a muoversi ognuno
%indipendentemente dall'altro.
%per calcolare il paramentro d'ordine posso riciclare il calcolo che ho
%fatto numericamente con il codice kuramototest tuttavia utilizzando una
%nuova funzione poichè ho bisno che dall'equzione differenziale emerga la
%dipendenza dal parametro d'ordine. In realtà ho bisogno soltanto di avere
%le soluzioni degli N oscillatori al variare del tempo perchè conosco una
%relazione del parametro d'ordine in funzione unicamente dei theta_i
%questo procedimento dovrei farlo due volte per i due coefficienti di
%accoppiamento.
tspan=[0,100];
delta_t=linspace(tspan(1),tspan(2),1500);
%delta_t=10;
theta_sol=[];
for i=1:length(K)
    kuraN = @(t,theta) (kura(t,theta,omega,K(i),N))';
    soluzione(i)=ode45(kuraN,tspan,iniz(:,i));
    theta_sol=[theta_sol; deval(delta_t,soluzione(i))];
end

%adesso farei un altro ciclo però per trovare il paramentro d'ordine
%quando richiamo delle componenti di matrici utilizzo le parentesi tonde
%per avere in ordine righe virgola colonne
%r=zeros(length(K),length(delta_t)); 
r=[];
for i=1:length(delta_t)
        rer1=(1/N)*sum(cos(theta_sol((1:N),i)));
        imr1=(1/N)*sum(sin(theta_sol((1:N),i)));
        rer2=(1/N)*sum(cos(theta_sol((N+1:2*N),i)));
        imr2=(1/N)*sum(sin(theta_sol((N+1:2*N),i)));
        r(1,i)=sqrt(rer1^2+imr1^2);
        r(2,i)=sqrt(rer2^2+imr2^2);
end
%il procedimento che sto qui utilizzando è di manipolare le righe della
%matrice soluzione theta sol in quanto per come è stata costruita un numero
%di righe pari al doppio degli oscillatori che ho supposto in partenza.
%Questo accade perchè il vettore k è lungo due e per ogni k desidero avere
%un set di soluzioni diverso in modo poi da contenere in r due righe quanti
%sono i k
plot(delta_t,r(1,:));
hold on;
plot(delta_t,r(2,:));
axis([0,tspan(2),0,1])
xlabel('t')
ylabel('r')
legend('k<k_c','k>k_c')
%attento agli indici i e j perchè il seno è dispari quindi invertendo gli
%indici su cui sommare nella fuznione si raggiungerebbero delle soluzioni 
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
%mi serve per ricavare da ode45 le soluzioni degli oscillatori nel modo più
%semplice possibile evitando di risolvere il sistema di equazioni
%differenziali apparentemente disaccoppiate ma con i parametri r e psi che
%dipendono dall'interazione di tutti gli oscillatori.

