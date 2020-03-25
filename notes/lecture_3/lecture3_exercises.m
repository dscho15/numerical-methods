%%  Numerical Methods Lecture 3

null([ 1 8 4 2 -4 ; 
       1 1 5 7 9  ;
       8 -5 1 -4 2;
       2 3 -1 3 -21;
       4 2 3 8 -8])
%% 

A = [ 1 0 1 ; 1 1 2 ; 0 1 1 ]
A_t = A.'

rref(A_t * A)

null(A_t * A)


(A_t * A)*20*[ 1.3 1 -1].'

%%

A = [ 1 1 2 ; 1 2 3 ; 1 2 3];  

c = [ 0 ; 0 ; 0]

for i = 1 : 1000

v = rand(3,1)

c = [ c A * v] 

end

for i = 1 : 100

v = 1*rand(1,1)*[ 1 1 -1].';

c = [ c v] 

v = 1*rand(1,1)*[1 1 1].';

c = [ c v] 

v = 1*rand(1,1)*[1 2 2].';

c = [ c v] 

v = 1*rand(1,1)*[2 3 3].';

c = [ c v] 

end

plot3(c(1,:), c(2,:), c(3,:),'*')


%%

A = [ 1 2 5 ; 2 4 10 ];

c = [ 0 0 0].'

basis = [ 1 2 5].';
u1 = [ 1 0 -1/5].';
u2 = [ 0 1 -2/5].';

for i = 1:100
   
    c = [c rand(1,1)*basis];
    c = [c rand(1,1)*u1];
    c = [c rand(1,1)*u2];
    
end


A * (u1)


plot3(c(1,:), c(2,:), c(3,:),'*')






