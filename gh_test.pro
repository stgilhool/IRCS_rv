function cb_hermite, order, x_in, sigma, norm

;n is order, x is the (x-x0) term in the Gaussian, and sigma is the
;sigma of the Gaussian, in the same units as x_in

n=order*1d0

x=x_in/sigma

h=0d

m=n/2l


for i=0l, m do h=h+(-1d0)^(i*1d0)*(2d0*x)^(n-2d0*i)/(factorial(i)*factorial(n-2*i))

output=norm*factorial(n)*h*exp(-x_in^2d0/(2d0*sigma^2d0))

return, output

end



function sg_hermite, order, x_in, sigma, norm

n=order*1d0

x=x_in/sigma

h=0d

m=n/2l

case order of
    0: h=1d0
    1: h=2*x
    2: h=4*x^2-2
    3: h=8*x^3-12*x
    4: h=16*x^4-48*x^2+12
    5: h=32*x^5-160*x^3+120*x
    6: h=64*x^6-480*x^4+720*x^2-120
    7: h=128*x^7-1344*x^5+3360*x^3-1680*x
    8: h=256*x^8-3584*x^6+13440*x^4-13440*x^2+1680
    9: h=512*x^9-9216*x^7+48384*x^5-80649*x^3+30240*x
    10: h=1024*x^10-23040*x^8+161280*x^6-403200*x^4+302400*x^2-30240
endcase

output=norm*h*exp(-x_in^2d0/(2d0*sigma^2d0))



return, output

end

function sg_hermite2, order, x_in, sigma, coeff, norm

n=order*1d0

x=x_in/sigma

h=0d

m=n/2l



if order eq 0 then h=1d0 else h=poly(x,coeff[*,order])

output=norm*h*exp(-x_in^2d0/(2d0*sigma^2d0))



return, output

end



function sg_hermite3, order, x_in, sigma, coeff

n=order*1d0

x=x_in/sigma

h=0d

m=n/2l



if order eq 0 then h=1d0 else h=poly(x,coeff[*,order])

output=h*exp(-x_in^2d0/(2d0*sigma^2d0))



return, output

end




pro gh_test

x=dindgen(51)-25
sig=5d0


coeff=dblarr(11, 11)
coeff[1,1]=2d0
coeff[2,2]=4d0
coeff[0,2]=-2d0
coeff[1,3]=-12d0
coeff[3,3]=8d0
coeff[0,4]=12d0
coeff[2,4]=-48d0
coeff[4,4]=16d0
coeff[1,5]=120d0
coeff[3,5]=-160d0
coeff[5,5]=32d0
coeff[0,6]=-120d0
coeff[2,6]=720d0
coeff[4,6]=-480d0
coeff[6,6]=64d0
coeff[1,7]=-1680d0
coeff[3,7]=3360d0
coeff[5,7]=-1344d0
coeff[7,7]=128d0
coeff[0,8]=1680d0
coeff[2,8]=-13440d0
coeff[4,8]=13440d0
coeff[6,8]=-3584d0
coeff[8,8]=256d0
coeff[1,9]=30240d0
coeff[3,9]=-80649d0
coeff[5,9]=48384d0
coeff[7,9]=-9216d0
coeff[9,9]=512d0
coeff[0,10]=-30240d0
coeff[2,10]=302400d0
coeff[4,10]=-403200d0
coeff[6,10]=161280d0
coeff[8,10]=-23040d0
coeff[10,10]=1024d0



for ord=0l, 10l do begin

norm_cb=(sig * 2^ord * factorial(ord) * sqrt(!pi))^(-0.5)
norm_sg=(sig^2 * 2^(ord+1) * factorial(ord) * (!pi))^(-0.5)
norm_sg2=(sig^2 * 2^(ord+1) * factorial(ord+1)* (!pi))^(-0.5)

;;;CB
    cb_tic=tic()
    
    lsf_cb=cb_hermite(ord, x, sig, norm_cb)
    

    cb_t=toc(cb_tic)
    
    plot, x, lsf_cb, title="Order : "+strtrim(ord,2)
    print, "CB Total : ", total(sqrt(lsf_cb^2))
    print, "CB Time : ", cb_t
    
;;;SG
    ;help, ord
    sg_tic=tic()
    lsf_sg=sg_hermite(ord, x, sig, norm_sg)
    sg_t=toc(sg_tic)
    
    oplot, x, lsf_sg, color=200
    print, "SG Total : ", total(sqrt(lsf_sg^2))
    print, "SG Time : ", sg_t

;;;SG2
    ;help, ord
    sg2_tic=tic()
    lsf_sg2=sg_hermite2(ord, x, sig, coeff, norm_sg2)
    sg2_t=toc(sg2_tic)
    
    oplot, x, lsf_sg2, color=99999
    print, "SG2 Total : ", total(sqrt(lsf_sg2^2))
    print, "SG2 Time : ", sg2_t
    

;;;;Gauss
    if ord eq 0 then begin
        gauss=(sig^2 * 2 * (!pi))^(-0.5)*exp(-1.*(x^2/(2*sig^2)))
        oplot, x, gauss, ps=6
    endif else begin
        comp=lsf_cb/(total(abs(lsf_cb)))
        print, "Comp Total : ", total(abs(comp))
        oplot, x, comp, ps=6
    endelse

;;;Final try
    lsff_tic=tic()
    lsff=sg_hermite3(ord, x, sig, coeff)
    lsffnorm=lsff/total(abs(lsff))
    lsff_t=toc(lsff_tic)
    oplot, x, lsffnorm, color=9999

    print, "final Total : ", total(sqrt(lsffnorm^2))
    print, "final Time : ", lsff_t
    

    stop
    
endfor

end
