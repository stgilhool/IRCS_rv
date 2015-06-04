function sg_hermite_coeff

coeff=dblarr(11, 11)
coeff[0,0]=1d0
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


return, coeff

end
