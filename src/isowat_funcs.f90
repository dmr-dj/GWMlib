!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2010 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module contient les variables commune pour les isotopes de 
!        l'eau. Module commun a priori aux parties océan et atmosphère
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 07 décembre 2010
!      Derniere modification : 23 juin 2021, 13 juillet 2023 (version 23), sept. 08th, 023
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       MODULE isowat_funcs

       use isowat_topdefs, only: ip, dblp=>dp
       use isowat_topdefs, only: nwisos, iwat18, iwat2h, iwat16, iwat17
       use isowat_defs,    only: M16, rsmow

       implicit none

       private :: nwisos, iwat18, iwat2h, iwat16, iwat17, M16, rsmow

       public :: rmois2nbmolesw, delta2moles, check_isowat_content, moles2delta, R2delta, delta2R, liq_vap_E, vap_liq_C

!      Définitions de compatibilité ...
        INTEGER(kind=ip), PARAMETER :: iso_noc = 1_ip, iso_nld = 3_ip, iso_nse = 2_ip


!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ ! dmr   Moved here the remnants of isoatm_mod ...
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        REAL(kind=dblp), DIMENSION(nlat, nlon, nwisos) :: ratio_oceanatm
!~        REAL(kind=dblp), DIMENSION(nlat, nlon, nwisos) :: ratio_evap_ocean, ratio_evap_land, ratio_evap_snow

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ ! dmr 2023-07-13  Restart flags, an odd place to have them, need to be updated at some point
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        INTEGER(kind=ip),PARAMETER :: isoatm_restart = WISOATM_RESTART
!~        INTEGER(kind=ip),PARAMETER :: isolbm_restart = WISOLND_RESTART


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! 
!       LIST OF THE FUNCTION WITHIN THIS LIBRARY

! 		function dexcess(d,o)             : computes the d-excess from the 2H ratio (d) and the 18O ratio (o), result ~ 1.E-3 (no per mil) 
!       function delta_invR(diso,n)       : computes the relative ratio using the delta of the given isotope
!       function deltaR(riso,n)           : computes the classical delta from the relative molar ratio (e.g. 18n/16n) for the given isotope      
!       delta2moles(rmois,deltaws)        : computes the number of moles for isotopes based on water content and conventionnal deltas
!       moles2delta(nb_moles)             : computes the conventional deltas (~ 1.E-3) given all isotopes molar content
!       evap_isoE(rmois, ratiosIso,tempK) : computes the molar content in the different water types for evaporating water (17->end)


!       function check_isowat_content(rwat, nmoleswiso) : check function for consistency between total water and isotopes
!       function REAL_EQUAL(realleft, realright) : function to test equality in the real sense, to an epsilon precision
!
! [NOTA] : REAL_EQUAL could be moved somewhere else, it is of more general interest!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Routine pour le calcul du deuterium-excess
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        pure elemental function dexcess(d,o) result (h2excess)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        d : ratio deuterium
!        o : ratio 18O
!
!       Variables de sortie : 
!        dexcess : le d-excess
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          IMPLICIT NONE

          REAL(kind=dblp)             :: h2excess
          REAL(kind=dblp), INTENT(IN) :: d,o

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          h2excess=(d/rsmow(iwat2h)-1.0_dblp)-8.0_dblp*(o/rsmow(iwat18)-1.0_dblp)

        end function dexcess

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Routine pour le calcul du ratio a partir du delta
!-----|--1--------2---------3---------4---------5---------6---------7-|
        pure elemental function delta2R(diso,n) result (riso_val)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        
!        diso : delta of the isotope considered
!        n : isotope number
!
!       Variables de sortie : 
!        riso_val: the relative molar ratio for the isotope considered
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       REAL(kind=dblp)             :: riso_val
       REAL(kind=dblp), INTENT(IN) :: diso
       integer(kind=ip), intent(in):: n

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       riso_val = rsmow(n) * (diso + 1._dblp) 

       end function delta2R
              
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Routine pour le calcul du delta a partir du ratio
!-----|--1--------2---------3---------4---------5---------6---------7-|
        pure elemental function R2delta(riso,n) result (delta_val)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        ratio : ratio in moisture
!        n : isotope number
!
!       Variables de sortie : 
!        delta : delta
!-----|--1--------2---------3---------4---------5---------6---------7-|
      
       IMPLICIT NONE

       REAL(kind=dblp)             :: delta_val
       REAL(kind=dblp), INTENT(IN) :: riso
       integer(kind=ip), intent(in):: n

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       delta_val = 0._dblp

       delta_val = (riso/rsmow(n) - 1._dblp)

       end function R2delta


      pure elemental function rmois2nbmolesw(rmois) result(nbmolesw)

      REAL(kind=dblp) :: nbmolesw
      REAL(kind=dblp), INTENT(IN) :: rmois

      ! dmr input rmoisg is in m3.m-2, output is in kmoles

      nbmolesw = ( rmois *1000.) / M16 ! M16 = molar mass of H216O 

      end function rmois2nbmolesw

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fonction pour calculer le contenu en moles isotopiques à partir du contenu en eau totale et du delta conventionnel
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       pure function delta2moles(rmois,deltaws) result (nb_moles)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        rmois : le contenu en eau en m3.m-2 (typiquement rmoisg)
!        deltaws = deltaw17, deltaw18, deltaw2h, les deltas conventionnels ( ~ 1.E-3 )
!
!       Variables de sortie : 
!        nb_moles : le tableau des contenus molaires pour 16->nwisos
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       IMPLICIT NONE

       REAL(kind=dblp), DIMENSION(iwat16:nwisos) :: nb_moles
       REAL(kind=dblp), DIMENSION(iwat17:nwisos), INTENT(in) :: deltaws
       REAL(kind=dblp), INTENT(IN) :: rmois

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: nb_molesw
       REAL(kind=dblp), dimension(iwat17:nwisos) :: iso_ratios
       integer(kind=ip) :: iso

        nb_molesw = rmois2nbmolesw(rmois)*1000._dblp ! output of rmois2nbmolesw is in kmoles
        
        do iso=iwat17,nwisos
          iso_ratios(iso) = delta2R(deltaws(iso),iso)
        enddo

        nb_moles(iwat16) = nb_molesw / (1._dblp +  SUM(iso_ratios))

        do iso=iwat17,nwisos
          nb_moles(iso) = iso_ratios(iso) * nb_moles(iwat16)

        enddo
        
        nb_moles(iwat16) = nb_moles(iwat16) / 1000._dblp ! iwat16 is in kmoles while the rest is in moles
        
       end function delta2moles


       pure function moles2delta(nb_moles) result (deltaws)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        nb_moles : le tableau des contenus molaires pour 16->nwisos
!       Variables de sortie : 
!        deltaws = deltaw17, deltaw18, deltaw2h, les deltas conventionnels ( ~ 1.E-3 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       
       IMPLICIT NONE

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), intent(in) :: nb_moles       
              
       ! OUT
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: deltaws

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       
       INTEGER(kind=ip) :: iso
       
       ! Conventional delta is ^in_w/^16n_w / rsmow
       do iso=iwat17,nwisos
         deltaws(iso) = (nb_moles(iso)/nb_moles(iwat16)) / (rsmow(iso) * 1000._dblp ) -1._dblp ! 1000 for kmoles -> moles in iwat16
       enddo              
       
       end function moles2delta
       

       pure function moles2R(nb_moles) result (ratiows)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        nb_moles : le tableau des contenus molaires pour 16->nwisos
!       Variables de sortie : 
!        deltaws = deltaw17, deltaw18, deltaw2h, les deltas conventionnels ( ~ 1.E-3 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       
       IMPLICIT NONE

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), intent(in) :: nb_moles       
              
       ! OUT
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: ratiows

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       
       INTEGER(kind=ip) :: iso
       
       ! Conventional ratio is ^in_w/^16n_w
       do iso=iwat17,nwisos
         ratiows(iso) = nb_moles(iso)/(nb_moles(iwat16) * 1000._dblp) ! 1000 for kmoles -> moles in iwat16
       enddo              
       
       end function moles2R
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Fonction pour calculer le fractionnement isotopique à l'évaporation
!       Utilises : l'eau liquide à évaporer
!                  les contenus isotopiques de l'eau liquide à évaporer
!                  la température (en KELVIN)
!
!                  renvoie le nombre de moles d'eau évaporée dans chaque phase isotopique 17->fin
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        


       pure function liq_vap_E(rmois_vap, molesisoliq, tempK, rmois_liq) result(nb_molesvap)
       
       USE isowat_alphas, only: alpha_lv

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), INTENT(in) :: molesisoliq
       REAL(kind=dblp), INTENT(IN) :: rmois_vap, tempK
       REAL(kind=dblp), INTENT(IN), OPTIONAL :: rmois_liq
       

       ! OUT
       REAL(kind=dblp), DIMENSION(iwat16:nwisos) :: nb_molesvap
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        
       
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: alphaeq, ratiosIso
       REAL(kind=dblp) :: coef1, coef2, nbmoleswatervap
       INTEGER(kind=ip):: iso
       
       nbmoleswatervap = rmois2nbmolesw(rmois_vap)*1000.0_dblp ! output of rmois2nbmolesw is in kmoles
       
       ratiosIso(:) = moles2R(molesisoliq)
       alphaeq(:) = alpha_lv(tempK)

       coef1 = 1._dblp+SUM(ratiosIso(:))
       coef2 = 1._dblp+SUM(ratiosIso(:)/alphaeq(:))

       if (PRESENT(rmois_liq)) then
         do iso=iwat17,nwisos
           nb_molesvap(iso) = molesisoliq(iso)/alphaeq(iso) * (rmois_vap/rmois_liq) * (coef1/coef2)
         enddo
       else
         do iso=iwat17,nwisos
           nb_molesvap(iso) = ratiosIso(iso)/alphaeq(iso) * nbmoleswatervap/coef2 
         enddo
       endif
       
       nb_molesvap(iwat16) = (nbmoleswatervap-sum(nb_molesvap(iwat17:nwisos)))/1000._dblp
       
       end function liq_vap_E

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!        Same as liq_vap_E but for vap_liq_C
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        


       pure function vap_liq_C(rmois_liq, molesisovap, tempK, rmois_vap) result(nb_molesliq)
       
       USE isowat_alphas, only: alpha_lv

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), INTENT(in) :: molesisovap
       REAL(kind=dblp), INTENT(IN) :: rmois_liq, tempK
       REAL(kind=dblp), INTENT(IN), OPTIONAL :: rmois_vap
       

       ! OUT
       REAL(kind=dblp), DIMENSION(iwat16:nwisos) :: nb_molesliq
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        
       
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: alphaeq, ratiosIso
       REAL(kind=dblp) :: coef1, coef2, nbmoleswaterliq
       INTEGER(kind=ip):: iso
       
       nbmoleswaterliq = rmois2nbmolesw(rmois_liq)*1000.0_dblp ! output of rmois2nbmolesw is in kmoles
       
       ratiosIso(:) = moles2R(molesisovap)
       alphaeq(:) = 1._dblp/alpha_lv(tempK)

       coef1 = 1._dblp+SUM(ratiosIso(:))
       coef2 = 1._dblp+SUM(ratiosIso(:)/alphaeq(:))

       if (PRESENT(rmois_vap)) then
         do iso=iwat17,nwisos
           nb_molesliq(iso) = molesisovap(iso)/alphaeq(iso) * (rmois_liq/rmois_vap) * (coef1/coef2)
         enddo
       else
         do iso=iwat17,nwisos
           nb_molesliq(iso) = ratiosIso(iso)/alphaeq(iso) * nbmoleswaterliq/coef2 
         enddo
       endif
       
       nb_molesliq(iwat16) = (nbmoleswaterliq-sum(nb_molesliq(iwat17:nwisos)))/1000._dblp
       
       end function vap_liq_C

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Fonction pour calculer le fractionnement isotopique à l'évaporation
!       Utilises : l'eau liquide à évaporer
!                  les contenus isotopiques de l'eau liquide à évaporer
!                  la température (en KELVIN)
!
!                  renvoie le nombre de moles d'eau évaporée dans chaque phase isotopique 17->fin
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        


       pure function sol_vap_S(rmois_vap, molesisoliq, tempK, rmois_liq, KF) result(nb_molesvap)
       
       ! in this routine, liq means solid ... !
       
       USE isowat_alphas, only: alpha_sv, alpha_sve

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), INTENT(in) :: molesisoliq
       REAL(kind=dblp), INTENT(IN) :: rmois_vap, tempK
       REAL(kind=dblp), INTENT(IN), OPTIONAL :: rmois_liq
       LOGICAL,         INTENT(IN), OPTIONAL :: KF
       

       ! OUT
       REAL(kind=dblp), DIMENSION(iwat16:nwisos) :: nb_molesvap
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        
       
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: alphaeq, ratiosIso
       REAL(kind=dblp) :: coef1, coef2, nbmoleswatervap
       INTEGER(kind=ip):: iso
       
       nbmoleswatervap = rmois2nbmolesw(rmois_vap)*1000.0_dblp ! output of rmois2nbmolesw is in kmoles
       
       ratiosIso(:) = moles2R(molesisoliq)
       if (PRESENT(KF)) then
         alphaeq(:) = alpha_sve(tempK)
       else
         alphaeq(:) = alpha_sv(tempK)
       endif 
       
       coef1 = 1._dblp+SUM(ratiosIso(:))
       coef2 = 1._dblp+SUM(ratiosIso(:)/alphaeq(:))

       if (PRESENT(rmois_liq)) then
         do iso=iwat17,nwisos
           nb_molesvap(iso) = molesisoliq(iso)/alphaeq(iso) * (rmois_vap/rmois_liq) * (coef1/coef2)
         enddo
       else
         do iso=iwat17,nwisos
           nb_molesvap(iso) = ratiosIso(iso)/alphaeq(iso) * nbmoleswatervap/coef2 
         enddo
       endif
       
       nb_molesvap(iwat16) = (nbmoleswatervap-sum(nb_molesvap(iwat17:nwisos)))/1000._dblp
       
       end function sol_vap_S

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!        Same as sol_vap_S but for vap_sol_C
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        


       pure function vap_sol_C(rmois_liq, molesisovap, tempK, rmois_vap, KF) result(nb_molesliq)

       ! in this routine, liq means solid ... !
       
       USE isowat_alphas, only: alpha_sv, alpha_sve

       REAL(kind=dblp), DIMENSION(iwat16:nwisos), INTENT(in) :: molesisovap
       REAL(kind=dblp), INTENT(IN) :: rmois_liq, tempK
       REAL(kind=dblp), INTENT(IN), OPTIONAL :: rmois_vap
       LOGICAL,         INTENT(IN), OPTIONAL :: KF

       ! OUT
       REAL(kind=dblp), DIMENSION(iwat16:nwisos) :: nb_molesliq
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|        
       
       REAL(kind=dblp), DIMENSION(iwat17:nwisos) :: alphaeq, ratiosIso
       REAL(kind=dblp) :: coef1, coef2, nbmoleswaterliq
       INTEGER(kind=ip):: iso
       
       nbmoleswaterliq = rmois2nbmolesw(rmois_liq)*1000.0_dblp ! output of rmois2nbmolesw is in kmoles
       
       ratiosIso(:) = moles2R(molesisovap)
       
       if (PRESENT(KF)) then       
         alphaeq(:) = 1._dblp/alpha_sve(tempK)
       else
         alphaeq(:) = 1._dblp/alpha_sv(tempK)
       endif
       
       coef1 = 1._dblp+SUM(ratiosIso(:))
       coef2 = 1._dblp+SUM(ratiosIso(:)/alphaeq(:))

       if (PRESENT(rmois_vap)) then
         do iso=iwat17,nwisos
           nb_molesliq(iso) = molesisovap(iso)/alphaeq(iso) * (rmois_liq/rmois_vap) * (coef1/coef2)
         enddo
       else
         do iso=iwat17,nwisos
           nb_molesliq(iso) = ratiosIso(iso)/alphaeq(iso) * nbmoleswaterliq/coef2 
         enddo
       endif
       
       nb_molesliq(iwat16) = (nbmoleswaterliq-sum(nb_molesliq(iwat17:nwisos)))/1000._dblp
       
       end function vap_sol_C

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fonction pour vérifier la coherence du contenu isotopique et de l'eau totale
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
       pure function check_isowat_content(rwat, nmoleswiso) result(OKISH)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        rwat : le contenu en eau en m3.m-2 (typiquement rmoisg)
!        nmoleswiso : le nombre de moles pour chaque isotope (iwat16 en kmoles pour conserver la précision) 
!
!       Variables de sortie : 
!        OKISH : un booléen vrai/faux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       IMPLICIT NONE
       
        REAL(kind=dblp), intent(in) :: rwat
        REAL(kind=dblp), dimension(nwisos), intent(in) :: nmoleswiso

        logical :: OKISH
        
        real(kind=dblp) :: nmolesright, nmolesleft
        
        nmolesleft = rmois2nbmolesw(rwat)

        nmolesright = nmoleswiso(iwat16) + SUM(nmoleswiso(iwat17:nwisos))/1000._dblp
        
        OKISH = REAL_EQUAL(nmolesleft, nmolesright)
       
       end function
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fonction pour calculer une égalité dans le monde réel avec une précision d'eps_dp
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
       pure elemental function REAL_EQUAL(realleft, realright) result(OKISH)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        realleft : le réel du membre de gauche
!        realright : le réel du membre de droite 
!
!       Variables de sortie : 
!        OKISH : un booléen vrai/faux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       
       use isowat_topdefs, only: eps_dp
       
       IMPLICIT NONE
       
       REAL(kind=dblp), intent(in) :: realleft, realright
       logical :: OKISH
        
        OKISH = .FALSE.
        OKISH = (ABS(realleft-realright).LE.eps_dp*1000._dblp)       
       
       end function REAL_EQUAL

       END MODULE isowat_funcs
