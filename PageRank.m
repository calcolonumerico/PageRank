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
%L'accuratezza dipende dal numero massimo di iterazioni e dalla tolleranza che sono settati 
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
%     20
%     19
%     19
%     20
%     20
%     .
%     .
% IN =
% 
%      0
%     18
%     16
%     16
%    .
%    .
%--------------------------------
% [U,G]=surfer('http://www.unina.it',50);
% [R, OUT, IN] = PageRank(G)
% 
% R =
% 
%    0.063631444211992
%    0.011005910489649
%    0.010721130453865
%    0.011005910489649
%    
% OUT =
% 
%     32
%      0
%     32
%      0
%      
% IN =
% 
%     34
%      3
%      2
%      3
% Autori:
%       Iodice Ivano
%       Vincenzo De Francesco

 if(nargin~=1)
     error("Inserire come parametro di input la matrice G")
 end

 %Controlli sulla matrice
    if(~ismatrix(G))
       error("Il primo input deve essere una matrice.")
    elseif(~issparse(G)||isempty(G))
       error("La matrice deve essere sparsa e non vuota.")
    elseif(size(G,1)~=size(G,2))
        error("La matrice deve essere quadrata.")
    elseif(~islogical(G))
      error("Gli elementi della matrice devono essere logici.")
    end

    p=0.85;
    n=size(G,1);
    e=ones(n,1);
   
    %Azzero la diagonale per evitare autoloop
    G(logical(eye(n)))=0;
   
    %Calcolo outdegree e indegree
    c = sum(G,1);
    OUT=full(c)';
    IN=full(sum(G,2));
    
    %Inizializzo z
    z = ((1-p)*(c~=0) + (c==0))/n; 
    
    %Inizializzo D
    D=sparse(1,n);
    D(c~=0)=1./c(c~=0);
    
    %Inizializzo Rank
    R=e/n;
    R0=zeros(n,1);
    
    niter=0;
    while(niter<200 && (norm(R-R0,Inf))>(10^-7)*norm(R0,Inf))
        R0=R;
        R=p*G.*D*R0+ones(n,1)*z*R0;
        niter=niter+1;
    end

end

