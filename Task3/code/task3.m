
% Calculate the first Model1 
load('resultsp1k6.mat');
load('PB12.mat','X1','X2');
x = vertcat(X1,X2);
 

[n D] = size(x); % number of observations (n) and dimension (D)
k = 6;           % number of components
  
clear Z;
   
for i=1:k
S1(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));
                 
end

S1 = S1./(2*p);
S1 = sum(S1,2);
    
 
%calculate the first Model2 
load('resultsp2k6.mat');
load('PB12.mat','X1','X2');
x = vertcat(X1,X2);
 
[n D] = size(x);     % number of observations (n) and dimension (D)
k = 6;               % number of components

 
% Calculate second model  
clear Z;
   
for i=1:k
S2(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));
                 
end
S2 = S2./(2*p);
S2 = sum(S2,2);
    
 
%Compare the different 
r_result = S1 > S2;
   
for i = 1:n
if(i<(n/2+1))
compare(i)=true;
else
compare(i)=false;
end
compare = compare';
end
Save_Tor = confusionmat(compare,r_result);
error = sum(compare~=r_result)/i;
   
       
