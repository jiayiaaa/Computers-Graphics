function Q=quaternion_mys(q,t)
% q=[2 1 1 3];t=2.2;
% q是单位四元数；
% t是单位四元数的幂；
q=dwh(q);

if  q(1)==1 
    Q=[1 0 0 0];
else
theta=acos(q(1));
m=sin(theta);
n=(1/m).*[q(2) q(3) q(4)];
A=cos(t.*theta);
B=sin(t.*theta).*n;
Q=cat(2,A,B);
end
end



%% 单位四元数的幂次运算结果应该也是个单位四元数