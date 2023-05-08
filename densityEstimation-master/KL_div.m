function div = KL_div(u1, S1, u2, S2)
% Compute the KL divergence between two Gaussians 
% KL(N(u1, S1)||N(u2, S2)) =
% 0.5*(Tr(S1^(-1)S2-I)-logdet(S2S1^(-1))+(u1-u2)'S1^(-1)(u1-u2)).
d = size(S1,1);
div = trace(S1^(-1)*S2 - eye(d,d));
div = div - log(det(S2*S1^(-1)));
div = div + (u1-u2)'*S1^(-1)*(u1-u2);
div = 0.5*div;