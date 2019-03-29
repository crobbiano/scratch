function [outs] = varyLambda(D,V,G,lambda)
K=(eye(size(D))+(1/lambda)*V'*V);
P1=(1/lambda^2)*D*V*G'*G*V'*D';
P2=(1/lambda^3)*D*V*K*V'*V*G'*G*V'*D';
P3=(1/lambda^3)*D*V*G'*G*V'*V*K*V'*D';
P4=(1/lambda^4)*D*V*K*V'*V*G'*G*K*V'*D';

outs = P1-P2-P3+P3;