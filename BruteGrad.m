% [X1 Y1 1;...
%  X2 Y2 1;...
%  X3 Y3 1];

Plane = [1 4 1;...
         2 3 1;...
         1 8 1];
Gradient(:,1) = Plane\[1;0;0];
Gradient(:,2) = Plane\[0;1;0]; 
Gradient(:,3) = Plane\[0;0;1];
Gradient(3,:) = 0; %Remove unwanted c coeff

x=1;
y=2;

n1(x) = Plane(1,1); %x1
n1(y) = Plane(1,2); %y1
n2(x) = Plane(2,1); %x2
n2(y) = Plane(2,2); %y2
n3(x) = Plane(3,1); %x3
n3(y) = Plane(3,2); %y3


if n2(y)~=n3(y)
    gn1(x) = -1/( ( (n2(x)-n3(x)) * (n1(y)-n3(y)) / (n2(y)-n3(y)) ) - (n1(x)-n3(x)) );
else
    gn1(x) = 0;
end
if n2(x)~=n3(x)
    gn1(y) = -1/( ( (n2(y)-n3(y)) * (n1(x)-n3(x)) / (n2(x)-n3(x)) ) - (n1(y)-n3(y)) );
else
    gn1(y) = 0;
end

if n1(y)~=n3(y)
    gn2(x) = -1/( ( (n1(x)-n3(x)) * (n2(y)-n3(y)) / (n1(y)-n3(y)) ) - (n2(x)-n3(x)) );
else
    gn2(x) = 0;
end
if n1(x)~=n3(x)
    gn2(y) = -1/( ( (n1(y)-n3(y)) * (n2(x)-n3(x)) / (n1(x)-n3(x)) ) - (n2(y)-n3(y)) );
else
    gn2(y) = 0;
end

if n1(y)~=n2(y)
    gn3(x) = -1/( ( (n1(x)-n2(x)) * (n3(y)-n2(y)) / (n1(y)-n2(y)) ) - (n3(x)-n2(x)) );
else
    gn3(x) = 0;
end
if n1(x)~=n2(x)
    gn3(y) = -1/( ( (n1(y)-n2(y)) * (n3(x)-n2(x)) / (n1(x)-n2(x)) ) - (n3(y)-n2(y)) );
else
    gn3(y) = 0;
end

% C style syntax
% if (n2[yy]!=n3[yy]){
%     gn1[xx] = -1/( ( (n2[xx]-n3[xx]) * (n1[yy]-n3[yy]) / (n2[yy]-n3[yy]) ) - (n1[xx]-n3[xx]) );
% }else{
%     gn1[xx] = 0;}
% if (n2[xx]!=n3[xx]){
%     gn1[yy] = -1/( ( (n2[yy]-n3[yy]) * (n1[xx]-n3[xx]) / (n2[xx]-n3[xx]) ) - (n1[yy]-n3[yy]) );
% }else{
%     gn1[yy] = 0;}
% 
% if (n1[yy]!=n3[yy]){
%     gn2[xx] = -1/( ( (n1[xx]-n3[xx]) * (n2[yy]-n3[yy]) / (n1[yy]-n3[yy]) ) - (n2[xx]-n3[xx]) );
% }else{
%     gn2[xx] = 0;}
% if (n1[xx]!=n3[xx]){
%     gn2[yy] = -1/( ( (n1[yy]-n3[yy]) * (n2[xx]-n3[xx]) / (n1[xx]-n3[xx]) ) - (n2[yy]-n3[yy]) );
% }else{
%     gn2[yy] = 0;}
% 
% if (n1[yy]!=n2[yy]){
%     gn3[xx] = -1/( ( (n1[xx]-n2[xx]) * (n3[yy]-n2[yy]) / (n1[yy]-n2[yy]) ) - (n3[xx]-n2[xx]) );
% }else{
%     gn3[xx] = 0;}
% if (n1[xx]!=n2[xx]){
%     gn3[yy] = -1/( ( (n1[yy]-n2[yy]) * (n3[xx]-n2[xx]) / (n1[xx]-n2[xx]) ) - (n3[yy]-n2[yy]) );
% }else{
%     gn3[yy] = 0;}