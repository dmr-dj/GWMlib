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
       LOGICAL, DIMENSION(iwat17:nwisos) :: test_result
            
       ! Testing fractionnation factors in delta mode
             
             
       delta_test = diniCste_m12pm
       ratio_test(iwat16:nwisos) = delta2R(delta_test(iwat16:nwisos),iwatns(:))

       WRITE(*,*) "Testing setup: "
       WRITE(*,*)  
       WRITE(*,*)  
       WRITE(*,*) "Initial deltas "
       WRITE(*,*) "============== "
       WRITE(*,1234) "d17O == ", delta_test(iwat17)
       WRITE(*,1234) "d18O == ", delta_test(iwat18)
       WRITE(*,1234) "d2H  == ", delta_test(iwat2h)
       WRITE(*,*)
       WRITE(*,*) "Fractionnation temperature for l/v "
       WRITE(*,*) "================================== "
       WRITE(*,1234) "tK_frac [K] = ", tK_frac
       WRITE(*,*)
       
!~        do temp = tK_zero_C,tK_zero_C+40.0d0
!~           alpha_test(iwat17:nwisos) = alpha_lv(temp)
!~           alpha_test2(iwat17:nwisos) = alpha_sv(temp-40.)
!~           write(*,*) temp, delta_test(iwat18)*1000._dblp,                                                          &
!~                                          R2delta(ratio_test(iwat18)/alpha_test(iwat18),iwat18)*1000._dblp&
!~                                                         ,alpha_test(iwat18), alpha_test2(iwat18)
!~        enddo     
              
              
       rmoisgtest = max_nb_rmoisg/2.
       moleswater = rmois2nbmolesw(rmoisgtest)
       

       WRITE(*,*) "Water content "
       WRITE(*,*) "============= "
       WRITE(*,1234) "in m^3.m^-2 ", rmoisgtest
       WRITE(*,1234) "in kmoles   ", moleswater
       WRITE(*,*)
              
       molesisowater(iwat16:nwisos) = delta2moles(rmoisgtest,delta_test(iwat17:nwisos))

       WRITE(*,*) "Isotopic content in moles "
       WRITE(*,*) "========================= "
       WRITE(*,1234) "16O == ", molesisowater(iwat16)*1000._dblp
       WRITE(*,1234) "17O == ", molesisowater(iwat17)
       WRITE(*,1234) "18O == ", molesisowater(iwat18)
       WRITE(*,1234) "2H  == ", molesisowater(iwat2h)
       WRITE(*,*)

       I_M_OK = check_isowat_content(rmoisgtest,molesisowater)

       WRITE(*,*) "Coherence check of isotopic content "
       WRITE(*,*) "=================================== "
       WRITE(*,2345) "Below tolerance? [T/F]", I_M_OK
       WRITE(*,*)


       test_result(:) = REAL_EQUAL(moles2delta(molesisowater(iwat16:nwisos)),delta_test(iwat17:nwisos))
       
       WRITE(*,*) "Absolute error on delta calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       
!~        write(*,*) "Delta this time", moles2delta(molesisowater(iwat16:nwisos))
!~        write(*,*) "Delta check", delta_test(iwat17:nwisos)

!~        liq_vap_E(rmoisgtest, ratio_test(iwat17:nwisos),tK_zero_C+20._dblp)
       write(*,*) "Equilibrium liquid vapor fractionnation"
       WRITE(*,*) "============== "


1234   FORMAT(4X, A, 4X, F12.3)       
2345   FORMAT(4X, A, 4X, L3)       
       end program main
