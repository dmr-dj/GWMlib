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
!
!      Auteur : Didier M. Roche 
!      Date   : 10 décembre 2010
!      Derniere modification : 14 décembre 2010, 23 juin 2021, Jul. 13h, 2023 (version 23), sept. 08th, 2023
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      module isowat_alphas

       use isowat_topdefs, only: ip, dblp=>dp, tK_zero_C
       use isowat_topdefs, only: nwisos, iwat18, iwat2h, iwat17, iwater
       use isowat_defs   , only: Diff

       implicit none

       private 
       
       public:: alpha_lv, alpha_sv, alpha_sve


          CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! 
!       LIST OF THE PUBLIC FUNCTIONS WITHIN THIS LIBRARY
!         
! 		function alpha_lv       : computes the relative fractionnation factor for liquid -- vapor, using the given temp in [K]
!       function alpha_sv       : id. for solid -- vapor
!       function alpha_lv       : id. for solid -- vapor with additional out of equilibrium fractionnation
!         
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       pure function alpha_lv(temp) result(alph_lvs)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement (KELVINS !!!)
!
!       Variables de sortie : 
!        alpha_lv : coefficients de fractionnement liquide-vapeur pour tous les isotopologues considérés
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp 
       REAL(kind=dblp), dimension(iwat17:nwisos) :: alph_lvs
       

       alph_lvs(iwat17) = alpha_lv17(temp)
       alph_lvs(iwat18) = alpha_lv18(temp)
       alph_lvs(iwat2h) = alpha_lv2h(temp)
       
       end function alpha_lv

       pure function alpha_sv(temp) result(alph_svs)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement (KELVINS !!!)
!
!       Variables de sortie : 
!        alpha_sv : coefficients de fractionnement solide-vapeur pour tous les isotopologues considérés
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp 
       REAL(kind=dblp), dimension(iwat17:nwisos) :: alph_svs
       

       alph_svs(iwat17) = alpha_sv17(temp)
       alph_svs(iwat18) = alpha_sv18(temp)
       alph_svs(iwat2h) = alpha_sv2h(temp)
       
       end function alpha_sv       

       pure function alpha_sve(temp) result(alph_sves)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement (KELVINS !!!)
!
!       Variables de sortie : 
!        alpha_sve : coefficients de fractionnement solide-vapeur pour tous les isotopologues considérés
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp 
       REAL(kind=dblp), dimension(iwat17:nwisos) :: alph_sves
       

       alph_sves(iwat17) = alpha_sve17(temp)
       alph_sves(iwat18) = alpha_sve18(temp)
       alph_sves(iwat2h) = alpha_sve2h(temp)
       
       end function alpha_sve


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! 
!       LIST OF THE PRIVATE FUNCTIONS WITHIN THIS LIBRARY

!       Fractionation function from the classical Majoube measurments and follow-ups

! 		function alpha_lv18       : computes the relative 18O fractionnation factor for liquid -- vapor, using the given temp in [K]
!       function alpha_lv2h(temp) : id. for 2H
!       function alpha_lv17(temp) : id. for 17O

!       function alpha_sv18(temp) : same as alpha_lv18 for vaporisation (solid -- vapor)
!       function alpha_sv2h(temp) : id. for 2H
!       function alpha_sv17(temp) : id. for 17O

!       function alpha_sve18(temp): same as alpha_sv18 but including the kinetic fractionnation of Jouzel
!       function alpha_sve17(temp): id. for 17O
!       function alpha_sve2h(temp): id. for 2H

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       pure elemental function alpha_lv18(temp) result(alph_lv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lv : coefficient de fractionnement liquide-vapeur for H218O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_lv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b, c

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fractionnement de Majoube, 1971b, Journal de Chimie Physique,
!         Volume 10, pp. 1423-1436
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
                 
        a = 1137.0d0 / temp**2
        b = 0.4156d0 / temp
        c = 2.0667d-3
        alph_lv = EXP(a-b-c)
    
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END FUNCTION alpha_lv18
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       pure elemental function alpha_lv2h(temp) result(alph_lv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lv : coefficient de fractionnement liquide-vapeur for H218O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_lv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b, c

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Fractionnement de Majoube, 1971b, Journal de Chimie Physique,
!         Volume 10, pp. 1423-1436
!-----|--1--------2---------3---------4---------5---------6---------7-|
                 
        a = 24844.0d0 / temp**2
        b = 76.248d0 / temp
        c = 56.612d-3
        alph_lv = EXP(a-b+c)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END FUNCTION alpha_lv2h
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       pure elemental function alpha_lv17(temp) result(alph_lv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lv : coefficient de fractionnement liquide-vapeur for H218O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_lv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b, c

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Fractionnement de Majoube, 1971b, Journal de Chimie Physique,
!         Volume 10, pp. 1423-1436
!       Et ajout de Barkan & Luz, 2005 et Barkan & Luz, 2007 pour 17O
!-----|--1--------2---------3---------4---------5---------6---------7-|

        a = 1137.0d0 / temp**2
        b = 0.4156d0 / temp
        c = 2.0667d-3
        alph_lv = (EXP(a-b-c))**0.529d0 ! Amaelle dixit : 0.529 pour lv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END FUNCTION alpha_lv17
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        pure elemental function alpha_sv18(temp)  result(alph_sv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alph_sv : coefficient de fractionnement solide-vapeur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL :: alph_sv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fractionnement de Majoube, 1971a, Journal de Chimie Physique,
!         Volume 68, pp. 625-636
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
                 
        a = 11.839d0 / temp
        b = 0.028244
        alph_sv = EXP(a-b)
    
      END FUNCTION alpha_sv18


       pure elemental function alpha_sv2h(temp)  result(alph_sv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alph_sv : coefficient de fractionnement solide-vapeur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL :: alph_sv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fractionnement de Majoube, 1971a, Journal de Chimie Physique,
!         Volume 68, pp. 625-636
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
                 
        a = 16288.0d0 / temp**2
        b = 0.0934
        alph_sv = EXP(a-b)

      END FUNCTION alpha_sv2h


        pure elemental function alpha_sv17(temp)  result(alph_sv)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alph_sv : coefficient de fractionnement solide-vapeur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), intent(in) :: temp ! attention, temp doit être en Kelvins
       REAL :: alph_sv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) :: a, b

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Fractionnement de Majoube, 1971a, Journal de Chimie Physique,
!         Volume 68, pp. 625-636
!       Et ajout de Barkan & Luz, 2005 et Barkan & Luz, 2007 pour 17O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        a = 11.839d0 / temp
        b = 0.028244
        alph_sv = (EXP(a-b))**0.528d0

      END FUNCTION alpha_sv17

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       pure elemental function alpha_sve18(temp) result(alph_sve)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lve : coefficient fractionnement solide-vapeur efficace
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), INTENT(IN) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_sve

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) a, b, alpha_eq, esse

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       REAL(kind=dblp), PARAMETER :: lambda = 0.004 ! entre 2. et 4.E-03 cf. Risi

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       On récupère d'abord la version à l'équilibre
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       alpha_eq = alpha_sv18(temp)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Le paramètre "S" est définit à partir de Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       esse = 1.0_dblp - lambda * (temp-tK_zero_C) ! here temp should be C not K

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Formulation du fractionnement cinétique lors de la formation 
!       des cristaux de glace, d'après Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       alph_sve = (alpha_eq * esse)/(1.0_dblp+alpha_eq*(esse-1.0_dblp)*Diff(iwater)/Diff(iwat18))

       END FUNCTION alpha_sve18

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!                                                                     |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       pure elemental function alpha_sve17(temp) result(alph_sve)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lve : coefficient fractionnement solide-vapeur efficace
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), INTENT(IN) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_sve

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) a, b, alpha_eq, esse

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       REAL(kind=dblp), PARAMETER :: lambda = 0.004 ! entre 2. et 4.E-03 cf. Risi

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       On récupère d'abord la version à l'équilibre
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       alpha_eq = alpha_sv17(temp)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Le paramètre "S" est définit à partir de Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       esse = 1.0_dblp - lambda * (temp-tK_zero_C) ! here temp should be C not K

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Formulation du fractionnement cinétique lors de la formation 
!       des cristaux de glace, d'après Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       alph_sve = (alpha_eq * esse)/(1.0_dblp+alpha_eq*(esse-1.0_dblp)*Diff(iwater)/Diff(iwat17))

       END FUNCTION alpha_sve17

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!                                                                     |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       pure elemental function alpha_sve2h(temp) result(alph_sve)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        temp : temperature de fractionnement
!
!       Variables de sortie : 
!        alpha_lve : coefficient fractionnement solide-vapeur efficace
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), INTENT(IN) :: temp ! attention, temp doit être en Kelvins
       REAL(kind=dblp) :: alph_sve

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp) a, b, alpha_eq, esse

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       REAL(kind=dblp), PARAMETER :: lambda = 0.004 ! entre 2. et 4.E-03 cf. Risi

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       On récupère d'abord la version à l'équilibre
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       alpha_eq = alpha_sv2h(temp)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Le paramètre "S" est définit à partir de Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       esse = 1.0_dblp - lambda * (temp-tK_zero_C) ! here temp should be C not K

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Formulation du fractionnement cinétique lors de la formation 
!       des cristaux de glace, d'après Jouzel & Merlivat, 1984
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       alph_sve = (alpha_eq * esse)/(1.0_dblp+alpha_eq*(esse-1.0_dblp)*Diff(iwater)/Diff(iwat2h))

       END FUNCTION alpha_sve2h

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   [TODO] rewrite the functions below as elemental and check the use in different realms                                      |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       FUNCTION alpha_ka(surf_type,iso_nb)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       Variables d'entree  : 
!~ !        surf_typ : le type de surface
!~ !        iso_nb : l'isotope considéré
!~ !
!~ !       Variables de sortie : 
!~ !        alpha_ka : le fractionnement cinétique de l'isotope considéré
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        USE iso_param_mod, ONLY : ieau, Diff, ieaud, ieau17, ieau18,     &
!~            iso_nld, iso_noc, iso_nse

!~        IMPLICIT NONE

!~        INTEGER :: surf_type, iso_nb
!~        REAL :: alpha_ka

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        REAL :: n_ka


!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Calcul sur les surfaces terrestres
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        IF (surf_type.EQ.iso_nld) THEN ! evaporation au dessus de la Terre
!~ !         n_ka => calcul à faire à partir de Mathieu & Bariac, 1996
!~          alpha_ka = ( Diff(ieau) / Diff(iso_nb) ) ** n_ka

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Calcul sur les surfaces océaniques
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        ELSE IF (surf_type.EQ.iso_noc) THEN

!~ !       Ici la détermination est à faire avec Amaëlle Landais
!~ !       Dans le cas de l'18O, 1+epsilon_ka = 6.0 pour mille
!~ !        => alpha_ka(ieau18) = 1 + 6.0D-3
!~ !       En attendant ...
!~          alpha_ka = 1.0
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Calcul sur les surfaces banquises
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        ELSE IF (surf_type.EQ.iso_nse) THEN
!~ !       Comme dans atmphys0.f, on dirait que evapn(noc) = evapn(nse)
!~ !       En attendant ...
!~          alpha_ka = 1.0

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Calcul sur les surfaces inconnues !!
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        ELSE
!~         WRITE(*,*) "Le type de surface demandé n'est pas codé !! "
!~         WRITE(*,*) "   >>> STOP <<<   "
!~         STOP
!~        ENDIF

!~       END FUNCTION alpha_ka

      END MODULE isowat_alphas
