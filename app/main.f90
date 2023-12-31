!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2023 Didier M. Roche (a.k.a. dmr)
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



       program main
       
       ! dmr testing program to be developped
       
       use isowat_topdefs, only: dblp=>dp, tK_zero_C, iwat16, iwat17, iwat18, nwisos, iwat2h, iwater, iwatns 
       use isowat_defs,    only: diniCste_m12pm,  rsmow, diniCste_m00pm,diniCste_m30pm
       use isowat_alphas,  only: alpha_lv, alpha_sv, alpha_sve
       use isowat_funcs,   only: rmois2nbmolesw, delta2moles, check_isowat_content, moles2delta, REAL_EQUAL, delta2R, R2delta &
                                , liq_vap_E, moles2R, vap_liq_C, vap_sol_C, sol_vap_S, nbmolesw2rmois, dexcessX
         
       implicit none
       
       REAL(dblp), DIMENSION(nwisos) :: delta_test, ratio_test, alpha_test1, alpha_test2, molesisowater, molesisowaterfinal   &
                                      , delta_fracR, moles_frac, delta_fracX, ratio2check                                     &
                                      , delta_fracX1, delta_fracX2, delta_fracX3, ratiolokaal
       
       REAL(dblp)                    :: tK_fraclv = tK_zero_C + 15._dblp, temp, moleswater,rmoisgtest                         &
                                      , tK_fracsv = tK_zero_C - 15._dblp
       
       REAL(dblp)                    :: min_nb_rmoisg = 0._dblp, max_nb_rmoisg = 5.2E-2_dblp ! in moles of water ...
       REAL(dblp)                    :: prop_evap = 4._dblp/5._dblp
       
       REAL(dblp)                    :: nmolesinvapor, nmolesinliquid
       REAL(dblp), DIMENSION(nwisos) :: vaporisocontent, liquidisocontent, delta_vapX, delta_liqX
       
       LOGICAL                       :: I_M_OK
       LOGICAL, DIMENSION(iwat17:nwisos) :: test_result
       
       INTEGER                       :: frac
       REAL(dblp)                    :: frac_vap
            
       ! Testing fractionnation factors in delta mode
             
             
       LOGICAL :: ACCURACY_WRITING = .FALSE.      
             
             
       delta_test = diniCste_m12pm
       ratio_test(iwat16:nwisos) = delta2R(delta_test(iwat16:nwisos),iwatns(:))

       WRITE(*,*) "Testing setup: "
       WRITE(*,*)  
       WRITE(*,*)  
       WRITE(*,*) "Initial deltas "
       WRITE(*,*) "============== "
       WRITE(*,1234) "d17O == ", delta_test(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_test(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_test(iwat2h)*1000._dblp
       WRITE(*,*)
       WRITE(*,*) "Fractionnation temperature for l/v "
       WRITE(*,*) "================================== "
       WRITE(*,1234) "tK_frac [K] = ", tK_fraclv
       WRITE(*,*)
       WRITE(*,*) "Fractionnation temperature for s/v "
       WRITE(*,*) "================================== "
       WRITE(*,1234) "tK_frac [K] = ", tK_fracsv
       WRITE(*,*)
       WRITE(*,1234) "prop. evap: ", prop_evap
       
!~        do temp = tK_zero_C,tK_zero_C+40.0d0
!~           alpha_test1(iwat17:nwisos) = alpha_lv(temp)
!~           alpha_test2(iwat17:nwisos) = alpha_sv(temp-40.)
!~           write(*,*) temp, delta_test(iwat18)*1000._dblp,                                                          &
!~                                          R2delta(ratio_test(iwat18)/alpha_test1(iwat18),iwat18)*1000._dblp&
!~                                                         ,alpha_test1(iwat18), alpha_test2(iwat18)
!~        enddo     
                            
       rmoisgtest = max_nb_rmoisg/2._dblp
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

       IF (ACCURACY_WRITING) THEN

       WRITE(*,*) "Coherence check of isotopic content "
       WRITE(*,*) "=================================== "
       WRITE(*,2345) "Below tolerance? [T/F]", I_M_OK
       WRITE(*,*)

       ENDIF
       
       test_result(:) = REAL_EQUAL(moles2delta(molesisowater(iwat16:nwisos)),delta_test(iwat17:nwisos))
       
       IF (ACCURACY_WRITING) THEN

       WRITE(*,*) "Absolute error on delta calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       
       ENDIF
       
       ratio2check(iwat17:nwisos) = moles2R(molesisowater(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(ratio2check(iwat17:nwisos),ratio_test(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN

       WRITE(*,*) "Absolute error on ratio calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)

       ENDIF
       
       WRITE(*,*) "Equilibrium liquid -> vapor fractionnation "
       WRITE(*,*) "========================================== "
       
       IF (ACCURACY_WRITING) THEN
       WRITE(*,*) "              VERSION 1                    "
       ENDIF

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = alpha_sv(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test1(iwat17:nwisos),iwatns(iwat17:nwisos))
              
       
       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (liquid implicit in iwat16) and rmoisg of vapor       
       moles_frac(iwat16:nwisos) = liq_vap_E(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fraclv)
       
       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF
       
       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
       
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = liq_vap_E(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fraclv,rmois_liq=rmoisgtest)

       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Evaporated deltas "
       WRITE(*,*) "================= "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)


       WRITE(*,*) "Equilibrium vapor -> liquid fractionnation "
       WRITE(*,*) "========================================== "
       IF (ACCURACY_WRITING) THEN
       WRITE(*,*) "              VERSION 1                 "
       ENDIF

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = 1._dblp/alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = 1._dblp/alpha_sv(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test1(iwat17:nwisos),iwatns(iwat17:nwisos))

       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (vapor implicit in iwat16) and rmoisg of liquid
       moles_frac(iwat16:nwisos) = vap_liq_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fraclv)

       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
              
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = vap_liq_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fraclv,rmois_vap=rmoisgtest)

       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Condensed deltas "
       WRITE(*,*) "================ "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)






! --- WITHOUT KINETIC FRACTIONNATION







       WRITE(*,*) "Equilibrium solid -> vapor fractionnation "
       WRITE(*,*) "========================================== "
       IF (ACCURACY_WRITING) THEN       
       WRITE(*,*) "              VERSION 1                    "
       ENDIF

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = alpha_sv(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test2(iwat17:nwisos),iwatns(iwat17:nwisos))
              
       
       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (liquid implicit in iwat16) and rmoisg of vapor       
       moles_frac(iwat16:nwisos) = sol_vap_S(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv)
       
       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
       
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = sol_vap_S(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv,rmois_liq=rmoisgtest)

       IF (ACCURACY_WRITING) THEN
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Evaporated deltas "
       WRITE(*,*) "================= "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)


       WRITE(*,*) "Equilibrium vapor -> solid fractionnation "
       WRITE(*,*) "========================================== "
       IF (ACCURACY_WRITING) THEN       
       WRITE(*,*) "              VERSION 1                 "
       ENDIF

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = 1._dblp/alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = 1._dblp/alpha_sv(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test2(iwat17:nwisos),iwatns(iwat17:nwisos))

       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (vapor implicit in iwat16) and rmoisg of liquid
       moles_frac(iwat16:nwisos) = vap_sol_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv)

       IF (ACCURACY_WRITING) THEN       
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF
       
       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN       
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
       
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = vap_sol_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv,rmois_vap=rmoisgtest)

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Condensed deltas "
       WRITE(*,*) "================ "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)


! --- WITH KINETIC FRACTIONNATION



       WRITE(*,*) "Equilibrium solid -> vapor fractionnation + KINETIC "
       WRITE(*,*) "=================================================== "
       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*) "                   VERSION 1                        "
       ENDIF
       
       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = alpha_sve(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test2(iwat17:nwisos),iwatns(iwat17:nwisos))
              
       
       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (liquid implicit in iwat16) and rmoisg of vapor       
       moles_frac(iwat16:nwisos) = sol_vap_S(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv,KF=.TRUE.)
       
       IF (ACCURACY_WRITING) THEN              
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF
       
       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
       
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = sol_vap_S(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv     &
                                            ,rmois_liq=rmoisgtest,KF=.TRUE.)
       IF (ACCURACY_WRITING) THEN              
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

 
       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN              
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Evaporated deltas "
       WRITE(*,*) "================= "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)


       WRITE(*,*) "Equilibrium vapor -> solid fractionnation + KINETIC"
       WRITE(*,*) "================================================== "
       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,*) "                   VERSION 1                       "
       ENDIF

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = 1._dblp/alpha_lv(tK_fraclv)
       alpha_test2(iwat17:nwisos) = 1._dblp/alpha_sve(tK_fracsv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test2(iwat17:nwisos),iwatns(iwat17:nwisos))

       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (vapor implicit in iwat16) and rmoisg of liquid
       moles_frac(iwat16:nwisos) = vap_sol_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv,KF=.TRUE.)

       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,*) "              VERSION 2                 "
       ENDIF
       
       ! VERSION 2, using moles isowater and rmoisg of vapor & liquid        
       moles_frac(iwat16:nwisos) = vap_sol_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fracsv         &
                                            ,rmois_vap=rmoisgtest,KF=.TRUE.)

       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,2345) " Below tolerance? [T/F]", check_isowat_content(rmoisgtest*prop_evap,moles_frac)
       ENDIF

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       IF (ACCURACY_WRITING) THEN                     
       WRITE(*,*)       
       WRITE(*,*) "Absolute error on deltaeq calculations below tolerance? [T/F]"
       WRITE(*,2345) "17O == ", test_result(iwat17)
       WRITE(*,2345) "18O == ", test_result(iwat18)
       WRITE(*,2345) "2H  == ", test_result(iwat2h)
       WRITE(*,*)
       ENDIF

       WRITE(*,*) "Condensed deltas "
       WRITE(*,*) "================ "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)









       WRITE(*,*) "Computation of final isotopic content "
       WRITE(*,*) "  -> example of liquid condensation   "
       WRITE(*,*) "     with d18O                        "
       WRITE(*,*) "=================================================== "
       WRITE(*,*)

       ! USING THE OLD DEFINITION / METHODOLOGY
       alpha_test1(iwat17:nwisos) = 1._dblp/alpha_lv(tK_fraclv)
       delta_fracR(iwat17:nwisos) = R2delta(ratio_test(iwat17:nwisos)/alpha_test1(iwat17:nwisos),iwatns(iwat17:nwisos))

       ! USING THE NEW DEFINITION / METHODOLOGY
       
       ! VERSION 1, using only moles isowater (vapor implicit in iwat16) and rmoisg of liquid
       moles_frac(iwat16:nwisos) = vap_liq_C(rmoisgtest*prop_evap, molesisowater(iwat16:nwisos),tK_fraclv)

       delta_fracX(iwat17:nwisos) = moles2delta(moles_frac(iwat16:nwisos))
       test_result(iwat17:nwisos) = REAL_EQUAL(delta_fracR(iwat17:nwisos),delta_fracX(iwat17:nwisos))

       WRITE(*,*) "Initial deltas "
       WRITE(*,*) "============== "
       WRITE(*,1234) "d17O == ", delta_test(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_test(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_test(iwat2h)*1000._dblp
       WRITE(*,*)


       WRITE(*,*) "Condensed deltas "
       WRITE(*,*) "================ "
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)
       
       WRITE(*,*) "Left in vapor "
       WRITE(*,*) "============= "
       
       molesisowaterfinal(iwat16:nwisos) = molesisowater(iwat16:nwisos) - moles_frac(iwat16:nwisos)
       
       delta_fracX(iwat17:nwisos) = moles2delta(molesisowaterfinal(iwat16:nwisos))
       
       WRITE(*,1234) "  %  == ", 1.0-prop_evap
       WRITE(*,1234) "d17O == ", delta_fracX(iwat17)*1000._dblp
       WRITE(*,1234) "d18O == ", delta_fracX(iwat18)*1000._dblp
       WRITE(*,1234) "d2H  == ", delta_fracX(iwat2h)*1000._dblp
       WRITE(*,*)
       
       WRITE(*,*) "Performing a Rayleigh distillation"
       WRITE(*,*) "================================= "   
       
       WRITE(*,*) "Starting point: "
       WRITE(*,*) "=============== "    
       
       
       nmolesinvapor = moleswater
       nmolesinliquid = 0._dblp
       
       vaporisocontent(iwater) = nmolesinvapor
       vaporisocontent(iwat16:nwisos)  = molesisowater(iwat16:nwisos)
       
       liquidisocontent(:) = 0._dblp
       
       ! delta_vapX, delta_liqX
       
       delta_vapX(iwat17:nwisos) = moles2delta(vaporisocontent(iwat16:nwisos))
       
       WRITE(*,1234) "d18Ov = ", delta_vapX(iwat18)*1000._dblp
       
       ratiolokaal(iwat17:nwisos) = moles2R(vaporisocontent(iwat16:nwisos))       
       
       do frac = 100, 1, -1
         
         frac_vap = real(frac)/100._dblp
         
         nmolesinvapor = moleswater*frac_vap            ! this is how much should left after this step
         nmolesinliquid = vaporisocontent(iwater) - nmolesinvapor ! this is how much we condensate        

         if (.not. REAL_EQUAL(nmolesinliquid,0._dblp)) then

         liquidisocontent(iwater) = nmolesinliquid
         liquidisocontent(iwat16:nwisos) = vap_liq_C(nbmolesw2rmois(nmolesinliquid), vaporisocontent(iwat16:nwisos) &
                                                   ,tK_fraclv) !,rmois_vap=nmolesinvapor)

         I_M_OK = check_isowat_content(liquidisocontent(iwater),liquidisocontent(:))
        
         delta_fracX1(iwat17:nwisos) = moles2delta(liquidisocontent(iwat16:nwisos))
         
         vaporisocontent(:) = vaporisocontent(:) - liquidisocontent(:)
                 
         delta_fracX2(iwat17:nwisos) = moles2delta(vaporisocontent(iwat16:nwisos))       
         
         
         
         
         WRITE(*,3456) "RayleigFrac", frac_vap, nmolesinvapor, nmolesinliquid, delta_fracX1(iwat18)*1000._dblp      &
                                    , delta_fracX2(iwat18)*1000._dblp                                               &
                                    , delta_fracX2(iwat2h)*1000._dblp                                               &
                                 , dexcessX(vaporisocontent(iwat16),vaporisocontent(iwat18),vaporisocontent(iwat2h))&
                                 , vaporisocontent(iwater)                                                          &
                                 , ((ratiolokaal(iwat18)*frac_vap**(alpha_lv(tK_fraclv)-1._dblp)) & 
                                                               /rsmow(iwat18)-1._dblp)*1000._dblp &                    
                                 , check_isowat_content(nbmolesw2rmois(vaporisocontent(iwater)),vaporisocontent(:)) &
                                 , check_isowat_content(nbmolesw2rmois(liquidisocontent(iwater)),liquidisocontent(:))
!~                                     , delta_fracX2(iwat2h)-8._dblp*delta_fracX2(iwat18)

         vaporisocontent(iwater) = nmolesinvapor

         endif


       enddo

       
       
       
1234   FORMAT(4X, A, 4X, F12.3)       
2345   FORMAT(4X, A, 1X, L3)   
3456   FORMAT(4X, A, 9F10.3, 2(2X, L3))      
       end program main
