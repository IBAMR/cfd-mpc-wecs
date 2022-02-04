function [A,B,F,C] = augmented_matrix(Phi,Gamma,Lambda)

m = size(Phi,2);

A = [Phi,          Gamma, Gamma;
     zeros(1,m),   1,     0;
     zeros(1,m),   0,     1];
B = [Lambda;
     1;
     0];
F  = [Lambda;
      0;
      1];

m = size(A,2); 
C = zeros(3,m);
C(1,1) = 1;
C(2,2) = 1;
C(3,m-1) = 1;