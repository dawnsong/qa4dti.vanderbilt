function [sortA sortB]=sortBbyA(A,B)

   %sorts B according to order in A
   % [sortA sortB]=sortBbyA(A,B)
   
   [sortA, ind]=sort(A);
   sortB=B(ind);
