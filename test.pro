;Just a worthless test of idl nuances
function testfunc, str

strtemp=str
strtemp.test1=str.test1+1

return, strtemp
end

pro test

t1=dindgen(10)
t2=dindgen(5)
parinfo={test1:t1, test2:t2}

a=testfunc(parinfo)

stop

end
