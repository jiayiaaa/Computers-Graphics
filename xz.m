% ����ά��v����λ��Ԫ��Q��ת�õ�vv
% v=[1 0 0];Q=[sqrt(2)/2 0 sqrt(2)/2 0];xz(v,Q)
% 
% ans =
% 
%     0.0000         0   -1.0000

function cc=xz(c,Q)
if length(c)==3
     c=cat(2,0,c);
 else
     c=c;
end
 Q=dwh(Q);
 Qc=quatmultiply(Q,c);
 QQ=quatconj(Q);
 cc=quatmultiply(Qc,QQ);
 cc=[cc(2) cc(3) cc(4)];
%  bs(Q);% ����γ�������ת��Ӧ�á�2
%  yv=[v(2) v(3) v(4)];
%  b=['��',num2str(yv),'��ת��',num2str(vv)];
%  disp(b);
 
end
 
