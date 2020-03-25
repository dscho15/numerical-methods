% Gramhams Method
clc; clear;

e = [ [ 4 -5 1 -4 3 ].' [ 1 1 5 7 8 ].' [ 2 8 4 2 1].']

e(:,1) = e(:,1)/norm(e(:,1));

for i = 2 : size(e,2)
    
    a = [0 0 0 0 0].';
    
    for j = 1 : i-1
        [ j i]
        a = a + dot(e(:,i),e(:,j)) * e(:,j)
        
    end
    
    e(:,i) = e(:,i) - a;
    e(:,i) = e(:,i)/norm(e(:,i));
    
end

dot(e(:,3),e(:,2))
