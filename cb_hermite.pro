function cb_hermite, n, x_in, sigma

;n is order, x is the (x-x0) term in the Gaussian, and sigma is the
;sigma of the Gaussian, in the same units as x_in

n=n*1d0

x=x_in/sigma

h=0d

m=n/2l

for i=0l, m do h=h+(-1d0)^(i*1d0)*(2d0*x)^(n-2d0*i)/(factorial(m)*factorial(n-2*i))

output=(1./sqrt(sigma))*(2d0^n*factorial(n)*Sqrt(!pi))^(-0.5d0)*factorial(n)*h*exp(-x_in^2d0/(2d0*sigma^2d0))



return, output



end
