clear, clc;
func = @(x)(x(2)-x(3))^2;
non1con=@mycon;

x0 = [10,10,10];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0,0];
ub = [21,21,21];

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

[x,fval] = fmincon(func,x0,A,b,Aeq,beq,lb,ub,non1con,options)

