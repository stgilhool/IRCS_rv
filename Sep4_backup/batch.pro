pro batch, cmdnum

case cmdnum of

    1: rvtest_display4, ft0='Apr07_lsf2_brute', dims=[11L,41L]

    2: rvtest_display4, expnum=[0], /bopt

    3: rvtest_display4, expnum=[0], /copt

    4: rvtest_display4, expnum=[0], ft0='Apr21_lsf3_brute', /bopt
    
    5: rvtest_display4, expnum=[0], ft0='Apr22_lsf3_wide', /bopt

    6: rvtest_display4, expnum=[0], ft0='Apr22_lsf3_more', /bopt

    7: rvtest_display4, expnum=[0], ft0='Apr23_lsf3_many', /bopt, dims=[61L,101L]

endcase

end
        
