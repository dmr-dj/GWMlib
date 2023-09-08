!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2021 Didier M. Roche (a.k.a. dmr)
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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      This is the base definition module for water isotopes, migrated from the LUDUS environnement
!      Now part of the more generic GWMlib
!
!      Auteur : Didier M. Roche 
!      Date   : August, 22nd, 2023 (version 23)
!      Latest modifications : Sept. 08th, 2023
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       MODULE isowat_defs_mod
       
       use isowat_topdefs_mod, only: ip, dblp=>dp
       use isowat_topdefs_mod, only: nwisos, iwat18, iwat2h, iwat16, iwat17

       implicit none

       private 

!  DEFINITIONS DU MODELE ISOTOPIQUE v2023
! 
!  1. Les variables prognostiques isotopiques contiennent les nombres de moles
!
!     C est a dire que nous avons :
!        indice 1 == ieau   | l'eau normale d'ECBilt, exprimee en [m3.m-2] = [m], exception car unite d'origine
!        indice 2 == ieau16 | le nombre de moles de la molecule H216O, exprimee en [kmols.m-2]. 
!                             N.B. le kmols est la pour assurer un ordre de grandeur de la variable proche de l'eau et des autres isotopes
!        indice 3 à 5       | le nombre de moles des molecules isotopiques de faible abondance: H217O, H218O, 1H2HO. En moles pour l'ordre de grandeur
!
!
!

!  2. Abondances relatives
!
!     Tous les calculs de fractionnement sont effectués en abondances molaires totales, notées iX, ou i est le numero de l'isotope et X le rapport molaire total. 
!       C est à dire que : iX = in / Somme_i in
!     L'avantage de cette notation est la conservation lors du melange. 
!
!
!     Pour assurer une compatibilite de sortie et d entrée avec ce qui est fait usuellement les routines de calcul des rapports d'abondance relative, les fonctions
!       sont toujours présentes ci-dessous mais sont annotées "R". Les nouvelles fonctions usuelles sont annotees "X" 
!


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Definition of the reference usual relative abundances
!      These are now compatibility references.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        REAL(KIND=dblp), PARAMETER :: r18smow = 2005.2D-6
        REAL(KIND=dblp), PARAMETER :: r17smow = 379.9D-6
        REAL(KIND=dblp), PARAMETER :: r2hsmow = 155.76D-6

!--- dmr  In practice in version 21, we will use 1.0 for any one of the isotopes
        REAL(KIND=dblp), PARAMETER :: rsmow(nwisos) = [1.0_dblp,1.0_dblp,r17smow,r18smow,r2hsmow]

!      Molar masses of the Mendeleiv table

        REAL(KIND=dblp), PARAMETER    :: MH = 1.007825032_dblp, M2H = 2.014101778_dblp                                     &
            , MO16 = 15.994914622_dblp, MO17 = 16.9991315_dblp, MO18 = 17.9991604_dblp 

!      Molar masses of the different types of water

        REAL(KIND=dblp), PARAMETER   :: M18 = 2._dblp * MH + MO18, M17 = 2. * MH + MO17, M16 = 2._dblp * MH + MO16         &
                                      , M19 = MH + M2H + MO16, M32 = MO16 + MO16, M33 = MO16 + MO17, M34 = MO16 + MO18

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      For kinetic fractionation, use diffusivity coefficient as defined
!      in reference document. Please note that the coefficient is in the
!      form [Di/D]^n, the coefficients given here are the Di/D
!      Given the formulation used (cf. ref. document), there is no such
!      coefficient for ieau or ieau16. Thus, they are set to 1.0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      Values are coming from Merlivat, 1978 & Barkan and Luz, 2007
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: Diff = [1.0d0, 1.0d0, 0.9855d0, 0.9723d0, 0.9755d0]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      For the special case of the ocean, we want the formulation to be
!      still in the type of [Di/D], though it is fixed in Merlivat and
!      Jouzel, 1979 that epsilon_k = alpha_k + 1 = 6 per mil for 18O.
!      Thus, n_ka_oc is fixed to ca. 0.214
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      However, it is also fixed in Merlivat and Jouzel, 1979 that
!      epsilon_k = alpha_k + 1 = 6 *0.88 per mil for D. 
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER, PRIVATE ::                                                                 &
         n_ka_oc = [1.0_dblp,1.0_dblp,LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18)),  LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18)), &
                                      LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18))]

!      From these considerations, it naturally folows that the kinetic
!      fractionation coefficients are, for the ocean:
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: alpha_diff_oc = Diff**n_ka_oc

!      Kinetic fractionnation coefficient for the land surface
!      evaporation

!      First guess
       REAL(kind=dblp), PARAMETER, PRIVATE           :: n_ka_lnd = 0.57_dblp

       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: alpha_diff_lnd = Diff**n_ka_lnd

       REAL(kind=dblp), PARAMETER                     :: fac_17Oexc = 0.528_dblp


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Initialisation constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr 2023-07-13   Deltas should be on the order of 1.E-3, not true everywhere at this date. 


       REAL(kind=dblp), PARAMETER :: d18iniCste12 = -12.E-3_dblp,                                               & 
          diniCste_m12pm(nwisos) = [0._dblp,0._dblp,(d18iniCste12+1._dblp)**fac_17Oexc-1.0_dblp, d18iniCste12,d18iniCste12*8._dblp]

       REAL(kind=dblp), PARAMETER :: d18iniCste30 = -30.E-3_dblp,                                               &
                    diniCste_m30pm(nwisos) = [0.d0,0.d0,(d18iniCste30/1000.+1.0d0)**0.528d0-1.0, d18iniCste30,d18iniCste30*8._dblp]

       END MODULE isowat_defs_mod
