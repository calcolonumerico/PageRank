%% Esecuzione dell'algoritmo PageRank
      load mathwork200.mat;
      [R, OUT, IN] = PageRank(G);
 %% Grafico a barre page rank
        figure (1)
        bar(R);
        title('Ranking delle pagine');
        xlabel('Pagine')
        ylabel('Percentuale di tempo di navigazione')
 %% Grafico struttura G e grafo associato
        figure(2)
        subplot(1,2,1)
        G=G.*~speye(size(G,1));
        spy(G);
        title('Matrice di connettività G');
        subplot(1,2,2)
        ff=digraph(G,'OmitSelfLoops');
        plot(ff,'NodeColor','y','LineWidth',1,'Marker','p','MarkerSize',2,'layout','force');
        title('Grafo associato a G');
        
        
%% Sottografo nodi rank > media rank
        G1=digraph(G,'OmitSelfLoops');
        G1.Nodes.Name=U;
        G1.Nodes.Rank=R;
        H=subgraph(G1,find(G1.Nodes.Rank > mean(R)));
        figure(3)
        plot(H,'NodeLabel',{},'NodeCData',H.Nodes.Rank,'Layout','force');
        title("Sottografo dei nodi con rank maggiore della media dei rank");
        colorbar
      
 %% Stampa dei primi 15 risultati per ranking con outdegree e indegree
        [R ,index]=sort(R,'descend');
        TOP15 = table(R(1:15),OUT(index(1:15)),IN(index(1:15)),'Rownames',U(index(1:15)),'VariableNames',{'Ranking' 'OUT' 'IN'});
        display(TOP15);

        