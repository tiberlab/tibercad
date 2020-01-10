function plot_unfolded_bandstructure(filename, Emin, Emax, dE, unfold)
% PLOT_UNFOLDED_BANDSTRUCTURE plot a band structure calculated 
%   in tibercad using the unfold_to flag
%
%

if (~exist('unfold', 'var'))
  unfold = 1;
end

d = importdata(filename);

k = d.data(:, 1);


% # of eigenvalues
Ne = (size(d.data, 2) - 1) / 2;



if (~exist('dE', 'var'))
  dE = 0.01;
end
  
E = Emin:dE:Emax;
    
A = zeros(length(E), length(k));

for ik=1:length(k)
  for ie=1:Ne
   Eki = d.data(ik, ie + 1);
   Wki = 1.0;
   if (unfold)
     Wki = d.data(ik, Ne + ie + 1);
   end
   
   A(:, ik) = A(:, ik) + Wki*lorentzian(E'-Eki);
   
  end
end

figure;
[X,Y] = meshgrid(k, E);
pcolor(X, Y, real(A))
%shading flat
shading('interp')
colormap('hot')
%m = colormap('gray');
%colormap(1-m);
%hold on
%plot(k, Epc, 'r');

%figure;
%hold on
%for ik=1:length(k)
%    scatter(k(ik)*ones(size(E)), E, 5, A(:, ik), 'filled');
%end
%m = colormap('gray');
%colormap(1-m);
end

function y = lorentzian(x)
g = 0.02;
y = 1.0/pi* g ./ (x.^2 + g^2);
end

function y = gauss(x)
g = 0.02;
y = 1.0/pi* exp(-x.^2 / g^2);
end
