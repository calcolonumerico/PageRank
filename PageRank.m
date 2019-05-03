function [R,OUT,IN] = PageRank(G)
%PageRank La funzione implementa l'algoritmo PageRank di Google per
%l'ordinamento dei siti web.
%
% Sintassi:
% [R,OUT,IN]=PageRank(G)
%
% Descrizione:
% [R,OUT,IN]=PageRank(G) calcola il ranking delle pagine della matrice di
% adiacenza G con i corrispondenti outdegree ed indegree.
%
% Parametri di ingresso:
%   G               = matrice sparsa di adiacenza relativa a d un
%   sottoinsieme del web.
%
% Parametri di uscita:
%   R               = vettore dei rank delle pagine.
%   OUT             = vettore dei corrispondenti outdegree.
%   IN              = vettore dei corrispondenti indegree.
%
% Diagnostica:
% Il programma si arresta mostrando un messaggio di errore se la matrice
% non è conforme alle specifiche del problema.
% 
% Accuratezza:
% L'accuratezza dipende dal numero massimo di iterazioni e dalla tolleranza che sono settati 
% rispettivamente a 200 e 10^-7.
%
% Algoritmo
% La funzione implementa l'algoritmo PageRank di Google.
%
% Esempi di utilizzo:
% load mathwork200.mat;
% [R, OUT, IN] = PageRank(G)
% 
% R =
% 
%    0.001759217077077
%    0.005536967249761
%    0.004787988545060
%    0.004798261903168
%    0.004798261903168
%    .
%    .
% OUT =
% 
%    (1,1)       20
%    (2,1)       19
%    (3,1)       19
%    (4,1)       20
%    (5,1)       20
%     .
%     .
% IN =
% 
%    (2,1)       18
%    (3,1)       16
%    (4,1)       16
%    (5,1)       16
%    .
%    .
%--------------------------------
% [U,G]=surfer('http://www.unina.it',50);
% [R, OUT, IN] = PageRank(G)
% 
% R =
% 
%    0.060476030042693 
%    0.009671544519183 
%    0.009671544519183 
%    0.009671544519183
%    .
%    .
% OUT =
% 
%   (1,1) 33
%   (3,1) 33
%   (7,1) 19
%   (8,1) 5
%   .
%   .
%      
% IN =
%
%   (1,1) 34
%   (2,1) 2
%   (3,1) 2
%   (4,1) 2
%   (5,1) 21
%   .
%   .
%
% Autori:
%       Iodice Ivano
%       Vincenzo De Francesco


%Controlli su input e output
 if(nargin~=1)
     error("Inserire come parametro di input la matrice G")
 end
if(nargout~=3)
    error("Inserire come parametri di output R, OUT, IN")
end
 %Controlli sulla matrice
 if(~ismatrix(G))
       error("Il primo input deve essere una matrice.")
 elseif(~issparse(G)||isempty(G))
       error("La matrice deve essere sparsa e non vuota.")
 elseif(size(G,1)~=size(G,2))
        error("La matrice deve essere quadrata.")
 elseif(size(G,1)<2)
         error("La matrice deve essere almeno 2x2")
  elseif(~islogical(G))
      error("Gli elementi della matrice devono essere logici.")
 end

    p=0.85;
    n=size(G,1);
    e=ones(n,1);
   
    %Azzero la diagonale per evitare autoloop
   G=G.*~speye(n);
   
    %Calcolo outdegree e indegree
    c = sum(G,1);
    OUT=c';
    IN=sum(G,2);
    
    %Inizializzo z
    z = ((1-p)*(c~=0) + (c==0))/n;
    
    %Inizializzo D
    D=sparse(1,n);
    D(c~=0)=1./c(c~=0);
    
    %Inizializzo Rank
    R0=e/n;
    R=p*G.*D*R0+e*(z*R0);
    niter=1;
    TOLX=10^-7;
    val=TOLX*norm(R,Inf);
    A=p*G.*D;
    while(niter<200 && (norm(R-R0,Inf))>TOLX*norm(R0,Inf))
        if(val>realmin)
                TOLX=val;
            else
                TOLX=realmin;
        end
        R0=R;
        R=A*R0+e*(z*R0);
        niter=niter+1;
    end
end
