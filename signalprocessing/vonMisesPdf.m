
function p = vonMisesPdf(theta,mu,k)

p = 1./(2*pi*besseli(0,k)) .* exp(k.*cos(theta - mu));