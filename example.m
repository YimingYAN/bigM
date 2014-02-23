clear all;
clc;


%% Exmaple 1 - check if max and <= work
%  max 2x1 + 3x2
%
%  s.t.  x1 + 4x2 <= 4
%       2x1 + 3x2 <= 6
%       2x1 +  x2 <= 4
%
%       x1>=0, x2>=0.

c = [ 2  3 ]';
A = [1 4; 2 3 ; 2 1];
b = [4 6 4]';
inq = [-1 -1 -1];

p = bigM(A,b,c,inq,'max'); p.solve; 

%% Example 2 - check if min and >= work
%  min -3x1 + x2 + x3
%
%  s.t.   x1 - 2x2 +  x3 <= 11
%       -4x1 +  x2 + 2x3 >= 3
%       -2x1 +        x3  = 1
%       
%       x1, x2, x3>=0.

c=[-3 1 1]';
b=[11 3 1]';
A=[1 -2 1;-4 1 2;-2 0 1];
inq=[-1 1 0];
p = bigM(A,b,c,inq,'min'); p.solve; 


%% Example 3 - check adding big M

c = [1 2 3]';
A = [0 1 2; 3 2 1];
b = [3 6]'; 
inq = [0 0 0];
p = bigM(A,b,c,inq,'min'); p.solve; 

%% Example 4 - check infeasible
%  max   x1 + 2x2
% s.t.   x1 +  x2 >= 4
%       3x1 + 2x2 <= 6
%        x1, x2 >=0
c = [1 2]';
A = [1 1 ; 3 2];
b = [4 6]';
inq = [1 -1];

p = bigM(A,b,c,inq,'max'); p.solve; 

%% Example 5 - check unbounded
%  min -x1-2x2
% s.t.    x1 + x2 >= 3
%         x1 - x2 >= 1
%       x1, x2 >= 0

c = [-1 -2]';
A = [1 1; 1 -1];
b = [3 ; 1];
inq = [1 1];

p = bigM(A,b,c,inq,'min'); p.solve; 

%% Example 6 - check infeasible with bigM

c = [ 1 1 1 1 ]';
A = [ 1 2 0 1; 0 0 1 3; 1 0 1 0 ];
b = [2 3 6 ]';
inq = [ -1 -1 1 ];
p = bigM(A,b,c,inq,'max'); p.solve; 

%% Example 7 - check unbounded with bigM
c = [ 1 2 3 4 ]';
A = [ 1 2 -1 2; 2 0 3 -1 ];
b = [ 4 6 ]';
inq = [ 1 0 ];
p = bigM(A,b,c,inq,'max'); p.solve; 




