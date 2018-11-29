% Initialise parameters
load('PB12.mat','X1','X2');
x = vertcat(X1,X2);
x = vertcat(X1,X2);

[n D] = size(x);

%Specifying the grid
[XX,YY] = meshgrid (x_axis,y_axis);

XX=XX';
XX=XX(:);


YY=YY';
YY=YY(:);


matrix=[XX,YY];
[n D] = size(matrix);

%Number of Cluster
k = 3;

maximum_a=max(x);
minimum_a=min(x);

x_axis = linspace(minimum_a(1),maximum_a(1));
y_axis = linspace(minimum_a(2),maximum_a(2));


%Model1

load('resultsp1.mat');

x = matrix;
[n D] = size(x);

%Number of cluster
k = 3;

clear Z;

for i=1:k

A1(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));

end

A1 = A1./(2*p);
A1 = sum(A1,2);

%Model2

load('resultsp2k3.mat');

x = matrix;
[n D] = size(x);

%Number of cluster
k = 3;

clear Z;

for i=1:k

B2(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));


end

B2 = B2./(2*p);
B2 = sum(B2,2);

I1 = A1 > B2;
m_last = vec2mat(I1,100);
imagesc(m_last);

