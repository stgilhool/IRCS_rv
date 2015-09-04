pro compare_models, model_1, model_2

readcol, 'model_info.dat', model_num, chi2, DoF, format="L,D,D"


;;;Print final results for a given file
index_1=where(model_num eq model_1)
index_1=index_1[0]
index_2=where(model_num eq model_2)
index_2=index_2[0]
chi2_1=chi2[index_1]
chi2_2=chi2[index_2]
DoF_1=DoF[index_1]
DoF_2=DoF[index_2]
stop
F_crit=((chi2_1-chi2_2)/(DoF_1-DoF_2))/(chi2_2/DoF_2)
prob=mpftest(F_crit, DoF_1-DoF_2, DoF_2)



print, "Replace model ", model_1," with model ", model_2
print, "F_ratio: ", F_crit
print, "Prob : ", prob

stop
end
