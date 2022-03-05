%Vogliamo scrivere un sistema del primo ordine della forma:
% (d/dt) theta=f(t,theta,omega,beta,n)
%dove theta è il vettore delle incognite, omega un vettore fissato mediante
%rand, e beta ed n sono due parametri scalari fissati una tantum all'inizio
%(n rappresenta il numero di incognite e pertanto di equazioni).
%Evidentemente l'informazione è contenuta nella funzione f, che viene
%codificata alla fine del codice mediante alcuni cicli; si noti che il
%tempo t va specificato come argomento anche se nella f non appare
%esplicitamente perché per usare ode45 ci serve una dipendenza temporale
%dal parametro t scalare libero. Alla fine del codice "produco" una
%funzione kura che restituisce il vettore f dei secondi membri delle varie
%equazioni in funzione di t, theta, omega, beta, n; nel codice dichiaro dei
%valori specifici per omega=rand..., beta=K con K prefissato, n=N con N già
%dichiarato, e costruisco una funzione secondaria kuraN che dipenda solo da
%t e theta e che numericamente assuma il valore "offerto" da kura valutando
%esplicitamente solo gli ultimi argomenti (quindi solo con t e theta
%indeterminati). Questo passaggio è importante perché ode45 richiede che la
%funzione che gli passiamo dipenda solo dal tempo e dal vettore delle
%incognite, ma non richiede necessariamente una dichiarazione (di kuraN in
%questo caso) a parte; difatti si può fare tutto in un colpo solo
%dichiarando in loco la funzione @(t,theta) kura(t,theta,omega,K,N) anziché
%semplicemente chiamando kuraN (ma in quest'ultimo caso evidentemente non 
%ho bisogno di specificare gli argomenti liberi e quelli già valutati).
%Entrambe le alternative equivalenti sono presenti nel codice (commentate o
%meno). Quello che segue è un utilizzo standard di ode45, deval per
%estrarre i valori numerici dalla struttura soluzione valutati ai tempi
%scelti, e un modo balordo per visualizzare il risultato.
%--------------------------------------------------------------------------
clear all
N=5;
K=2;
%k non deve variare in un intervallo ben preciso poichè rappresenta il
%coupling degli oscillatori più o meno intenso

omega=randn(N,1);
iniz=rand(N,1)*2*pi;

%dipende quale dei due metodi sto utilizzando per valutare numericamente le
%soluzioni del sistema di equazioni differenziali avrò il vettore della
%soluzione come riga oppure colonna rispettivamente con tante colonne o
%righe quanti sono gli oscillatori del sistema
%omega=randn(1,N); 

kuraN = @(t,theta) (kura(t,theta,omega,K,N))';
tspan=[0,10];

soluzione=ode45(kuraN, tspan, iniz);
%[T, theta_sol] = ode45(kuraN, tspan , omega);
delta_t=linspace(tspan(1),tspan(2),200);
theta_sol=deval(delta_t,soluzione);

%theta_sol=deval(T,soluzione);
%  for i=1:N
%       plot(delta_t,theta_sol(i,:));%riga i= theta_i(t)
%       hold on;
%  end
%questo primo plot è interessante per vedere come tra loro sono
%sincronizzate le fasi degli oscillatori. Cioè capiamo in quanti
%oscillatori si discostano dal gruppo al variare di k. Questo conto
%andrebbe confrontato con il valore di k critico che si trova
%analiticamente tramite l'approssimazione del continuo di kuramoto in base
%alle distribuzioni delle frequenze

theta = 0:0.01:2*pi;
%delta t dovrebbe essere il numero di volte che ho valutato la soluzione di
%ode come anche si vede, oppure T che valuta lui stesso automaticamente
figura = figure;
filename = 'test.gif';
for k=1:length(delta_t) %length(T)
    p1=[0,0];
    p2=[(1/N)*sum(cos((theta_sol(:,k)))),(1/N)*sum(sin((theta_sol(:,k))))];
    df=p2-p1;
    quiver(p1(1),p1(2),df(1),df(2),0)
    %è necessario per rapresentare una freccia il comando quiver il quale
    %utilizza due punti che gli fornisco in input la loro differenza e
    %infine un fattore di scala che è posto uguale a zero come richiesto
    %dal tutorial su mathworks
    hold on;
    plot(cos(theta),sin(theta),cos(theta_sol(:,k)),sin(theta_sol(:,k)),'b*');
    %calcolo il parametro d'ordine e plotto sulla circonferenza
    %goniometrica
    %plot(cos(theta),sin(theta),cos(theta_sol(k,:)),sin(theta_sol(k,:)),'*');
    axis equal
    axis([-1.5 1.5 -1.1 1.1])
    pause(0.1);
    drawnow limitrate;
    hold off;
    hold off;
    % Capture the plot as an image 
      frame = getframe(figura); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if k == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
%drawnow;
    %ATTENZIONE C'E UNO STRETTO LEGAME TRA PAUSE E TSPAN CIOE IL TEMPO NEL
    %QUALE IL COMANDO ODE INTEGRA LE EQUAZIONI E PAUSE. PRATICAMENTE IL
    %PROGRAMMA IMPIEGA UN CERTO TEMPO PER ESSERE ESEGUITO E QUINDI FARE I
    %PLOT. PAUSE INDICA CHE OGNI TOT SECONDI IL PROGRAMMA SI FERMA E MI FA
    %VEDERE IL PLOT IN QUEL MOMENTO COSì FACENDO IL PROGRAMMA è RALLENTATO.
    %SE IO INVECE DI USARE 1 SECONDO USASSI UN MILLISECONDO IL PROGRAMMA
    %VERREBBE,CREDO, ESEGUITO PIU VELOCEMENTE (COME NORMALMENTE DOVREBBE
    %ESSERE) VEDENDO COSì IL PLOT MUOVERSI PIU RAPIDAMENTE. PER TSPAN
    %ELEVATI NELLO STESSO TEMPO DI ESECUZIONE DEL PROGRAMMA O QUASI IL PLOT
    %DEVE PROCEDERE IN UN INTERVALLO MAGGIORE QUINDI SEMBRA SEMPRE PIù
    %VELOCE. PERTANTO DOVRO ALL' AUMENTARE DI TSPAN INCREMENTARE ANCHE IL
    %TEMPO DI PAUSE CAUSANDO PERO' UN' EVOLUZIONE PIU' DISCONTINUA DEL
    %SISTEMA.
     
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