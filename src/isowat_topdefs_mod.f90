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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      This is the top definition module for water isotopes, designed as a top level compatibility module for the library
!      Now part of the more generic GWMlib
!
!      Auteur : Didier M. Roche 
!      Date   : Sept. 08th, 2023 (version 23)
!      Latest modifications : 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       MODULE isowat_topdefs
       

       implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Selected double precision as taken from:
!   Metcalf, M., J. Reid, and M. Cohen (2004). Fortran 95/2003 Explained. Oxford University Press.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, parameter  :: sp = 4      ! kind(1.0) does not work if compilation is forced to double as in e.g. -fdefault-real-8
       integer, parameter  :: dp = selected_real_kind(2*precision(1.0_sp))
       integer, parameter  :: qp = selected_real_kind(2*precision(1.0_dp))

       integer, parameter  :: sip = 2, ip = 4, dip = 8

       integer, parameter  :: str_len =256
       integer, parameter  :: uuid_size = 36

       INTEGER(ip), parameter      :: nwisos = 5
       INTEGER(ip), parameter      :: iwater = 1
       INTEGER(ip), parameter      :: iwat16 = 2
       INTEGER(ip), parameter      :: iwat17 = 3
       INTEGER(ip), parameter      :: iwat18 = 4
       INTEGER(ip), parameter      :: iwat2h = 5       


! dmr  Temperature of 0°C in K                   (in K)
       real(dp), parameter :: tK_zero_C  = 273.15_dp

! afq  Tiny numbers
       real(sp), parameter :: eps_sp = epsilon(eps_sp)
       real(dp), parameter :: eps_dp = epsilon(eps_dp)

       END MODULE isowat_topdefs
