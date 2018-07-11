function y = voigt(x, sigma, gamma)

z = (x + 1i*gamma)/(sqrt(2) * sigma);
y = real(fadf(z))/(sqrt(2*pi) * sigma);