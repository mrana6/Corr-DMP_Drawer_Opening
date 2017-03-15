function h = plotGMM(Mu, Sigma, color, display_mode, dispSigma);
%
% This function plots a representation of the components (means and 
% covariance amtrices) of a Gaussian Mixture Model (GMM) or a
% Gaussian Mixture Regression (GMR).
%
% Author:	Sylvain Calinon, 2010
%			http://programming-by-demonstration.org/SylvainCalinon
%
% Inputs -----------------------------------------------------------------
%   o Mu:           D x K array representing the centers of the K GMM components.
%   o Sigma:        D x D x K array representing the covariance matrices of the 
%                   K GMM components.
%   o color:        Color used for the representation
%   o display_mode: Display mode (1 is used for a GMM, 2 is used for a GMR
%                   with a 2D representation and 3 is used for a GMR with a
%                   1D representation).

nbData = size(Mu,2);
lightcolor = color + [0.4,0.4,0.4];
lightcolor(find(lightcolor>1.0)) = 1.0;

if display_mode==1
  nbDrawingSeg = 20;
  t = linspace(-pi, pi, nbDrawingSeg)';
  for j=1:nbData
    stdev = sqrtm((dispSigma'*dispSigma).*Sigma(:,:,j));
    X = [cos(t) sin(t)] * real(stdev) + repmat(Mu(:,j)',nbDrawingSeg,1);
    h((j-1)*2+1)=patch(X(:,1), X(:,2), lightcolor, 'lineWidth', 2, 'EdgeColor', color);
    h((j-1)*2+2)=plot(Mu(1,j), Mu(2,j), 'x', 'lineWidth', 2, 'color', color);
    %[V,D]=eig(Sigma(:,:,j));
    %text(Mu(1,j), Mu(2,j), num2str(det(Sigma(:,:,j))));
    %text(Mu(1,j), Mu(2,j), [num2str(D(1,1)) ' , ' num2str(D(1,1))]);
    %prop=max([D(1,1) D(2,2)])/min([D(1,1) D(2,2)]);
    %text(Mu(1,j), Mu(2,j), num2str(prop));
  end
elseif display_mode==2
  nbDrawingSeg = 40;
  t = linspace(-pi, pi, nbDrawingSeg)';
  for j=1:nbData
    stdev = sqrtm((dispSigma'*dispSigma).*Sigma(:,:,j));
    X = [cos(t) sin(t)] * real(stdev) + repmat(Mu(:,j)',nbDrawingSeg,1);
    patch(X(:,1), X(:,2), lightcolor, 'LineStyle', 'none');
  end
  plot(Mu(1,:), Mu(2,:), '-', 'lineWidth', 1, 'color', color);
elseif display_mode==3
  for j=1:nbData
    ymax(j) = Mu(2,j) + sqrtm(dispSigma.*Sigma(1,1,j));
    ymin(j) = Mu(2,j) - sqrtm(dispSigma.*Sigma(1,1,j));
  end
  patch([Mu(1,1:end) Mu(1,end:-1:1)], [ymax(1:end) ymin(end:-1:1)], lightcolor, 'LineStyle', 'none');
  plot(Mu(1,:), Mu(2,:), '-', 'lineWidth', 3, 'color', color); 
elseif display_mode==4
  for j=1:nbData
    ymax(j) = Mu(2,j) + sqrtm(dispSigma.*Sigma(1,1,j));
    ymin(j) = Mu(2,j) - sqrtm(dispSigma.*Sigma(1,1,j));
  end
  verylightcolor = lightcolor + [0.2,0.2,0.2];
  verylightcolor(find(verylightcolor>1.0)) = 1.0;
  %'facealpha',0.5
  patch([Mu(1,1:end) Mu(1,end:-1:1)], [ymax(1:end) ymin(end:-1:1)], verylightcolor, 'LineStyle', 'none');
  plot(Mu(1,:), Mu(2,:), '--', 'lineWidth', 2, 'color', lightcolor); 
end





