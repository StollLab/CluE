function pass = test_findClusters()
   pass = true;
   side = 2;
   Nuclei = test_findClusters_generateNuc(side);
   cluster2 = findClusters(Nuclei,Nuclei.Index,2,1);
   % 1
   if norm(size(cluster2)-[12,2]) ~= 0
       pass = false;
       fprintf('\nfindClusters() --> failed at test 1.\n');
       return;
   end
   % 2
   cluster2 = findClusters(Nuclei,Nuclei.Index,2,1.5);
   if norm(size(cluster2)-[24,2]) ~= 0
       pass = false;
       fprintf('\nfindClusters() --> failed at test 2.\n');
       return;
   end
   % 3
   cluster2 = findClusters(Nuclei,Nuclei.Index,2,2);
   if norm(size(cluster2)-[28,2]) ~= 0
       pass = false;
       fprintf('\nfindClusters() --> failed at test 3.\n');
       return;
   end
   clear cluster2;
   % 4
   cluster3 = findClusters(Nuclei,Nuclei.Index,3,1);
   if norm(size(cluster3)-[24,3]) ~= 0
       cluster3
       size(cluster3)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 4.\n');
       return;
   end
   % 5
   cluster3 = findClusters(Nuclei,Nuclei.Index,3,1.5);
   if norm(size(cluster3)-[56,3]) ~= 0
       cluster3
       size(cluster3)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 5.\n');
       return;
   end
   % 6
   cluster3 = findClusters(Nuclei,Nuclei.Index,3,2);
   if norm(size(cluster3)-[56,3]) ~= 0
       size(cluster3)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 6.\n');
       return;
   end
   clear cluster3;
   % 6.1
   cluster4 = findClusters(Nuclei,Nuclei.Index,4,2);
   if norm(size(cluster4)-[70,4]) ~= 0
       size(cluster4)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 6.1.\n');
       return;
   end
   % 6.2
   cluster4 = findClusters(Nuclei,Nuclei.Index,4,1);
   if norm(size(cluster4)-[30,4]) ~= 0
       cluster4
       size(cluster4)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 6.2.\n');
       return;
   end
   clear cluster4;
   % 7: 3 by 3 by 3 nuclei
   side = 3;
   Nuclei = test_findClusters_generateNuc(side);
   cluster2 = findClusters(Nuclei,Nuclei.Index,2,1);
   if norm(size(cluster2)-[54,2]) ~= 0
       cluster2
       size(cluster2)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 7.\n');
       return;
   end
   clear cluster2;
   % 8
   cluster3 = findClusters(Nuclei,Nuclei.Index,3,1);
   if norm(size(cluster3)-[171,3]) ~= 0
       cluster3
       size(cluster3)
       pass = false;
       fprintf('\nfindClusters() --> failed at test 4.\n');
       return;
   end
   % Takes about a minute to run if test 9 is included.
   if false  
     clear cluster3;
     % 9: 10 by 10 by 10
     side = 10;
     Nuclei = test_findClusters_generateNuc(side);
     cluster4 = findClusters(Nuclei,Nuclei.Index,4,1);
     q = size(cluster4,1)/150000; % |clusters4| should be approximately 1000*6*5*5 = 150000, and is 47940, from code.
     log_q = log(q)/log(10);
     if  ~(log_q>0) && ~(log_q<1)
       q*150000
       pass = false;
       fprintf('\nfindClusters() --> failed at test 4.\n');
       return;
     end
   end
   fprintf('\nfindClusters() --> pass.\n');
end

function Nuclei = test_findClusters_generateNuc(side)
   numNuc=side^3;
   x = side;
   y = x;
   z = x;
   for ii = numNuc:-1:1
       Nuclei.Index(ii) = ii;
       z = mod(z-1, side);
       if z==(side-1)
           y = mod(y-1, side);
           if y==(side-1)
               x = mod(x-1, side);
           end
       end
       Nuclei.Coordinates(ii,:) = [x,y,z];
   end
end