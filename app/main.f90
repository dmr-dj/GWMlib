       program main
       
       ! dmr testing program to be developped
       
       use isowat_topdefs, only: dblp=>dp, tK_zero_C, iwat16, iwat17, iwat18, nwisos, iwat2h, iwater, iwatns 
       use isowat_defs,    only: diniCste_m12pm,  rsmow, diniCste_m00pm,diniCste_m30pm
       use isowat_alphas,  only: alpha_lv, alpha_sv
       use isowat_funcs,   only: rmois2nbmolesw, delta2moles, check_isowat_content, moles2delta, REAL_EQUAL, delta2R, R2delta &
                                , liq_vap_E
         
       implicit none
       
       REAL(dblp), DIMENSION(nwisos) :: delta_test, ratio_test, alpha_test, alpha_test2, molesisowater
       REAL(dblp)                    :: tK_frac = tK_zero_C + 15._dblp, temp, moleswater,rmoisgtest
       
       REAL(dblp)                    :: min_nb_rmoisg = 0._dblp, max_nb_rmoisg = 5.2E-2_dblp ! in moles of water ... 
       
       LOGICAL                       :: I_M_OK
     
     
       ! Testing fractionnation factors in delta mode
             
       delta_test = diniCste_m30pm
       ratio_test(iwat16:nwisos) = delta2R(delta_test(iwat16:nwisos),iwatns(:))
                    
       do temp = tK_zero_C,tK_zero_C+40.0d0
          alpha_test(iwat17:nwisos) = alpha_lv(temp)
          alpha_test2(iwat17:nwisos) = alpha_sv(temp-40.)
          write(*,*) temp, delta_test(iwat18)*1000._dblp,                                                          &
                                         R2delta(ratio_test(iwat18)/alpha_test(iwat18),iwat18)*1000._dblp&
                                                        ,alpha_test(iwat18), alpha_test2(iwat18)
       enddo     
              
       rmoisgtest = max_nb_rmoisg/2.
       moleswater = rmois2nbmolesw(rmoisgtest)
       
              
       molesisowater(iwat16:nwisos) = delta2moles(rmoisgtest,delta_test(iwat17:nwisos))

       write(*,*) "indexes :", iwater, iwat16, iwat17, iwat18, iwat2h
       write(*,*) "molesisowater", molesisowater
       
       I_M_OK = check_isowat_content(rmoisgtest,molesisowater)
       
       write(*,*) "Am I OKish ???", I_M_OK       
       write(*,*) "Delta this time", moles2delta(molesisowater(iwat16:nwisos))
       write(*,*) "Delta check", delta_test(iwat17:nwisos)
       write(*,*) "Relative error on delta", REAL_EQUAL(moles2delta(molesisowater(iwat16:nwisos)),delta_test(iwat17:nwisos))
       
       write(*,*) "Equilibrium liquid vapor fractionnation", liq_vap_E(rmoisgtest, ratio_test(iwat17:nwisos),tK_zero_C+20._dblp)
       
       end program main
