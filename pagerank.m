function [R,OUT,IN] = pagerank(G)
%pagerank La funzione implementa l'algoritmo PageRank di Google per
%l'ordinamento dei siti web
%   Detailed explanation goes here

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

    n=size(G,1);
    %Azzero la diagonale per evitare autoloop
    G=G-spdiags(spdiags(G,0),0,n,n);
    %Calcolo outdegree e indegree
    D=sum(G);
    OUT=full(D)';
    IN=full(sum(G,2));
    
    %Inizializzo z
    z=ones(1,n)*((1-0.85)/n);
    %Pongo z = 1/n per i dangling nodes
    z(D==0)=1/n;
    
    %Inizializzo D
    D(D~=0)=1./D(D~=0);
    D=diag(D);
    
    %Inizializzo Rank
    R=1/n*ones(n,1);
    R0=zeros(n,1);
    
    niter=0;
    while(niter<200 && (norm(R-R0,Inf))>(10^-7)*norm(R0,Inf))
        R0=R;
        R=0.85*G*D*R0+ones(n,1)*z*R0;
        niter=niter+1;
    end

end

