clear 
clf
%quello che si vuole fare in questo codice è rappresentare il parametro
%d'ordine in funzione della costante di accoppiamento. Significa che al
%posto di memorizzare le soluzioni questa volta posso utilizzarle subito
%per caclolare il parametro d'ordine per ogni k differente.
N=200; 
K=linspace(0,8);
%il coefficiente di accoppiamento andrà al massimo da zero a otto per un
%totale di 100 punti per rappresentare alla fine il grafico r,k in quanto i
%tempi tspan sono proprio 100
omega=tan(pi*(rand(N,1)-0.5)); 
%questa volta le frequenze sono distribuite come una lorentziana
iniz=zeros(N,1);
%le condizioni iniziali sono poste tutte uguali a zero per semplicità in
%quanto stiamo seguendo come anche suppone kuramoto nel caso stazionario al
%limite del continuo (caso che stiamo analizzando noi) le fasi devono
%essere stazionarie e nulla vieta di porsi in un sistema di riferimento
%corotante dove il parametro d'ordine nel TEMPO deve restare costante
%quindi indipendentemente da come partono gli oscillatori la loro coerenza
%deve restare la stessa; naturalmente cambierà al variare invece della
%costante di accoppiamento che è il secondo parametro, oltre il tempo
%rispetto al quale in questa situazione è stazionario, da cui dipende r.
tspan=[0,100];
delta_t=linspace(tspan(1),tspan(2));
%adesso procedo all'integrazione senza memorizzare i theta sol ma facendo
%dipendere il numero di gruppi di soluzioni per gli oscillatori da k.
r=[];
for i=1:length(K)
    kuraN = @(t,theta) (kura(t,theta,omega,K(i),N))';
    soluzione(i)=ode45(kuraN,tspan,iniz);
    theta_sol=deval(delta_t,soluzione(i));   
    %così ho ogni k una matrice per le soluzioni di cento righe cioè cento 
    %oscillatori. Ogni k dovrò trovare il parametro d'ordine corrispondente
    for j=1:length(delta_t)
     rer=(1/N)*sum(cos(theta_sol(:,j)));
     imr=(1/N)*sum(sin(theta_sol(:,j)));
     %rer e imr sono scalari essendo emersi da una somma 
     r(i,j)=sqrt(rer^2+imr^2);
    end
end
%vengono fuori 10000 righe eventualmente se memorizzassimo theta sol,questo
%significa che ogni cento righe è un valore di k ben preciso 
%voglio un vettore contenente le medie dei parametri d'ordine che ho
%trovato. Voglio un vettore con tante colonne quanti sono i parametri
%d'ordine quindi posso fare un ciclo che si muova sulle colonne e per
%ognuna di queste mi scriva il valor medio tra tutte le colonne di r (il
%ciclo sarà sulle righe invece di r)
%faccio tutto tramite il comando mean! questo comando, insieme ad errorbar
%restituisce, se applicato ad una matrice mxn una matrice riga con n
%colonne nelle quali sono contenute le medie eseguite per ogni n 
errorbar(K,mean(r'),std(r'),'.k');
hold on;
km = @(x) sqrt(1-2/x); 
fplot(km,[2,8]) 
axis([0,8,0,1]) 
xlabel('K') 
ylabel('r')
legend('numerico','teorico')
%dove teorico si intende l'evoluzione al variare di K del parametro
%d'ordine quando la ditribuzione delle frequenze è lorentziana unimodale.



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

