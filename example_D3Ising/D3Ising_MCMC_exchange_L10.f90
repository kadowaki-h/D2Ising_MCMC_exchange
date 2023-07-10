! D3Ising_MCMC_exchange_L10.f90  2023/7/4
!
! gfortran -fopenmp -O -o D3Ising_MCMC_exchange_L10.ex mt19937_64_OMP.f90 D3Ising_MCMC_exchange_L10.f90
! ifort -qopenmp -O -o D3Ising_MCMC_exchange_L10.ex mt19937_64_OMP.f90 D3Ising_MCMC_exchange_L10.f90
! export OMP_NUM_THREADS=16
! ./D3Ising_MCMC_exchange_L10.ex
!
!**********************************************************************
! Copyright (c) 2023, Hiroaki Kadowaki [email: Kadowaki.Hiroaki.G at gmail.com (remove spaces)]
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, 
! this list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation and/or 
! other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors may 
! be used to endorse or promote products derived from this software without specific 
! prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
! POSSIBILITY OF SUCH DAMAGE.
! 
!    (See https://opensource.org/licenses/BSD-3-Clause )
!
!**********************************************************************
!
! A fortran90 program of a MCMC simulation with the replica exchange algorithm
! which calculates <E(T)>, <C(T)>, <|M|(T)> etc of a 2d Ising model.
!-----
!  Markov chain Monte Carlo method
!  Metropolis algorithm
!  single spin flip
!  exchange Monte Carlo
!  Ising model
!  d=3 LxLxL cubic lattice
!  E = - Jnn * sum_{x = 1,L; y = 1,L; z = 1,L} 
!                [ spin(x,y,z)*spin(x+1,y,z) + spin(x,y,z)*spin(x,y+1,z) + spin(x,y,z)*spin(x,y,z+1) ]
!     spin(x,y,z) = 1 or -1
!     spin(L+1, y , z ) = spin(1,y,z)
!     spin( x ,L+1, z ) = spin(x,1,z)
!     spin( x , y ,L+1) = spin(x,y,1)
!  M = sum_{x=1,L; y=1,L; z=1,L} spin(x,y,z)
!  temperature (k_B * T) = temperature_list( 1, 2, ... , N_temperature_list )
!
!-----
MODULE exchange_MC
  use,intrinsic :: iso_fortran_env
  IMPLICIT NONE
  private
  public :: Jnn, L, N_spin, N_mcs_discard, N_mcs_average
  public :: R_num, R_T, R_spin, R_E, R_sum_observables
  public :: R_perm, R_accept_EXMC
  public :: R_mti, R_mt
  public :: sum_observable_data, output_table
!
  real(REAL64), save :: Jnn
  integer, save :: L, N_spin, N_mcs_discard, N_mcs_average
  integer, save :: R_num                                      !! number of replicas [temperatures (magnetic fields ... )]
  real(REAL64), save, allocatable :: R_T( : )                 !! temperatures R_T(j) > R_T(j+1) j = 1, 2, ...
  integer, save, allocatable :: R_spin( :,:,:,: )             !! spin of a replica R_spin(:,:,:, j ) j = 1, 2, ... , R_num
  integer, save, allocatable :: R_perm( : )                   !! R_T(j) corresponds to R_spin(:,:,:, R_perm(j) )
  real(REAL64), save, allocatable :: R_E( : )                 !! energy of a replica R_spin(:,:,:, j )
  real(REAL64), save, allocatable :: R_sum_observables( :,: ) !! summation E, E**2, M, M**2 etc R_sum_observables(:,j) ( at T = R_T(j) )
  real(REAL64), save, allocatable :: R_accept_EXMC( : )
!
  integer, save, allocatable :: R_mti( : )                    !! states of MT pseudo-random number generator
  integer(int64), save, allocatable :: R_mt( :,: )
!
  contains
!
!-----
SUBROUTINE sum_observable_data( j, E, M, accept )
  use,intrinsic :: iso_fortran_env
!!!!!  use exchange_MC , only : N_spin, R_sum_observables
!  summation of observables
  IMPLICIT NONE
  integer, intent(in) :: j
  real(REAL64), intent(in) :: E, M, accept
  real(REAL64) :: magnetization
  R_sum_observables(1, j) = R_sum_observables(1, j) + accept / DBLE(N_spin)  !! sum accept/N_spin
  R_sum_observables(2, j) = R_sum_observables(2, j) + E                      !! sum E
  R_sum_observables(3, j) = R_sum_observables(3, j) + E*E                    !! sum E**2
  magnetization = M / DBLE(N_spin)
  R_sum_observables(4, j) = R_sum_observables(4, j) + magnetization          !! sum M/N_spin
  R_sum_observables(5, j) = R_sum_observables(5, j) + abs(magnetization)     !! sum |M|/N_spin
  R_sum_observables(6, j) = R_sum_observables(6, j) + magnetization**2       !! sum (M/N_spin)**2
  R_sum_observables(7, j) = R_sum_observables(7, j) + magnetization**4       !! sum (M/N_spin)**4
END SUBROUTINE sum_observable_data
!
!-----
SUBROUTINE output_table( nf )
  use,intrinsic :: iso_fortran_env
!!!!!  use exchange_MC , only : N_spin, R_num, N_mcs_average, R_T, R_sum_observables, R_accept_EXMC
  IMPLICIT NONE
  integer, intent(in) :: nf
  real(REAL64) :: temperature, acceptance, energy_av, energy_sq_av, specific_heat
  real(REAL64) :: m_av, m_abs_av, m_sq_av, m_sqsq_av, chi, abs_chi, binder, accept_EXMC
  integer :: j
!
  do j = 1, R_num
    temperature   = R_T(j)
    acceptance    = R_sum_observables(1, j) / DBLE(N_mcs_average)
    energy_av     = R_sum_observables(2, j) / DBLE(N_mcs_average)  !! <E>
    energy_sq_av  = R_sum_observables(3, j) / DBLE(N_mcs_average)  !! <E**2>
    specific_heat = (energy_sq_av - energy_av * energy_av) / ((temperature**2) * DBLE(N_spin) )
!                                                                  !! <C/N_spin> = (1/T)**2 * (<E**2>-<E>**2)/N_spin
    energy_av     = energy_av / DBLE(N_spin)                       !! <E/N_spin>
    m_av          = R_sum_observables(4, j) / DBLE(N_mcs_average)  !! <M/N_spin>
    m_abs_av      = R_sum_observables(5, j) / DBLE(N_mcs_average)  !! <|M|/N_spin>
    m_sq_av       = R_sum_observables(6, j) / DBLE(N_mcs_average)  !! <(M/N_spin)**2>
    m_sqsq_av     = R_sum_observables(7, j) / DBLE(N_mcs_average)  !! <(M/N_spin)**4>
    chi     = DBLE(N_spin) * ( m_sq_av - (m_av**2)     ) / temperature !! (1/T) * ( <M**2> - <M>**2 )/ Nspin
    abs_chi = DBLE(N_spin) * ( m_sq_av - (m_abs_av**2) ) / temperature !! (1/T) * ( <M**2> - <|M|>**2 )/ Nspin
    binder  = m_sqsq_av / ( m_sq_av**2 )                               !! <M**4>/<M**2>**2
    accept_EXMC = R_accept_EXMC( j ) / ( 0.5d0 * DBLE(N_mcs_average) )
    write(nf,'(20G15.7)') &
    & temperature, energy_av, specific_heat, chi, abs_chi, m_av, m_abs_av, m_sq_av, binder, acceptance, accept_EXMC
  end do
END SUBROUTINE output_table
!
END MODULE exchange_MC
!
!-----
PROGRAM main
!
!$ use omp_lib
  use,intrinsic :: iso_fortran_env
  use exchange_MC
  use mt19937_64 , only : init_genrand64, genrand64_real3
  use mt19937_64 , only : nn_mt64 => nn, mti_mt64 => mti, mt_mt64 => mt
  IMPLICIT NONE
!
  integer, allocatable :: spin( :,:,: )
  integer :: nf, iter, j, N_temperature_list, N_observable, x, y, z
  real(REAL64) :: accept, E, M, temperature
  real(REAL64), allocatable :: temperature_list( : )
  character(128) :: file_table
  integer(int64) :: seed_64bit
!
! initial settings
    seed_64bit = 21478288371993_int64 + 202306251056_int64
!!!!!    seed_64bit = seed_64bit + 202306291046_int64       !! any 64-bit integer of the seed
    L = 10                                       !! LxLxL cubic lattice L = 10  16  32  64 ...
    Jnn = 1.0d0                                  !! exchange constant  Jnn > 0 ferromagnetic 
    file_table = "D3Ising_MCMC_exchange_L10.tab" !! output file, "D3Ising_MCMC_exchange_L16.tab" ...
    N_temperature_list = 26  ;  allocate( temperature_list( N_temperature_list ) )  !! number of temperatures
    temperature_list( 1:N_temperature_list ) = & !! temperature_list(j) > temperature_list(j+1) j = 1, 2, ...
      & (/                               5.5d0,    5.4d0, 5.3d0, 5.2d0, 5.1d0, 5.0d0 &
        &  , 4.9d0, 4.8d0, 4.7d0, 4.6d0, 4.5d0,    4.4d0, 4.3d0, 4.2d0, 4.1d0, 4.0d0 &
        &  , 3.9d0, 3.8d0, 3.7d0, 3.6d0, 3.5d0,    3.4d0, 3.3d0, 3.2d0, 3.1d0, 3.0d0 /)
    N_mcs_discard = 2000000 !! 500000            !! MC steps (per spin) to be discarded (for thermalization)
    N_mcs_average = 2000000 !! 500000            !! MC steps (per spin) for averaging
    N_spin = L**3                                !! number of spins
!
! open file_table and write several lines
    open( newunit = nf , file = TRIM(file_table) )
    write(nf,'(A,A)')     '  file_table         = ', TRIM(file_table)
    write(nf,'(A,5G15.7)')'  Jnn                = ', Jnn
    write(nf,'(A,5I10)')  '  L, N_spin          = ', L, N_spin
    write(nf,'(A,5I10)')  '  N_temperature_list = ', N_temperature_list
    write(nf,'(A)')       '  temperature_list   = '
    write(nf,'(10X,5G15.7)') temperature_list(1:N_temperature_list)
    write(nf,'(A,5I20)')  '  seed_64bit         = ', seed_64bit
    write(nf,'(A,5I10)')  '  N_mcs_discard, N_mcs_average = ', N_mcs_discard, N_mcs_average
    write(nf,'(A,A)') "  temperature  <E/N_spin>  <C/N_spin>  <chi/N_spin>  <abs_chi/N_spin>" &
    & , "  <M/N_spin>  <|M|/N_spin>  <(M/N_spin)**2>  Binder_<M**4>/<M**2>**2  acceptance  acceptance_exchange"
!
! allocate arrays
    allocate( spin( L,L,L ) )
    R_num = N_temperature_list
    allocate ( R_T( R_num ) )
    allocate ( R_spin( L,L,L, R_num ) )
    allocate ( R_perm( R_num ) )
    allocate ( R_E( R_num ) )
    allocate ( R_accept_EXMC( R_num ) )
    N_observable = 16  !! number of observables
    allocate ( R_sum_observables( N_observable, R_num ) )
    allocate ( R_mti( R_num + 1 ) )
    allocate ( R_mt( nn_mt64, R_num + 1 ) )
    R_T( 1:R_num ) = temperature_list( 1:N_temperature_list )  !! set temperatures R_T(j) = temperature_list(j)
    R_perm( 1:R_num ) = (/ ( j, j = 1, R_num ) /)              !! initialize R_perm(j) = j
!
!! initialize pseudo-random generator mt19937_64
!! set random spin configuration
!$OMP parallel do default(none) &
!$OMP private( j, x, y, z ) &
!$OMP shared( R_num, seed_64bit, L, R_spin, R_mti, R_mt ) schedule(dynamic,1)
  do j = 1, R_num
    call init_genrand64( ( seed_64bit + (79*j+1) ) )         !! initialize pseudo-random generator mt19937_64
    DO z = 1, L
      DO y = 1, L
        DO x = 1, L
          R_spin(x,y,z, j) = -1
          IF( genrand64_real3() < 0.5D0 ) R_spin(x,y,z, j) = 1  !! random spin configuration (replica j)
        END DO
      END DO
    END DO
    R_mti(j) = mti_mt64  ;  R_mt( : , j) = mt_mt64( : )        !! keep the mt19937_64 state (replica j)
  end do
  j = R_num + 1
    call init_genrand64( ( seed_64bit + (79*j+1) ) )           !! initialize pseudo-random generator mt19937_64
    R_mti(j) = mti_mt64  ;  R_mt( : , j) = mt_mt64( : )        !! keep the mt19937_64 state (for exchange MC)
!
! exchange MC, thermalization
  DO iter = 1, N_mcs_discard                                !! MC steps to be discarded (for thermalization)
!$OMP parallel do default(none) &
!$OMP private( j, temperature, spin, E, accept ) &
!$OMP shared( L, R_num, R_T, R_spin, R_E, R_perm, R_mti, R_mt ) schedule(dynamic,1)
    do j = 1, R_num
      temperature = R_T(j)
      mti_mt64 = R_mti(j)  ;  mt_mt64( : ) = R_mt( : , j)   !! restore the mt19937_64 state (replica j)
!
      spin(1:L,1:L,1:L) = R_spin(1:L,1:L,1:L, R_perm(j) )   !! restore the state of the replica
      accept = 0.0d0
        CALL one_MC_step( temperature, spin, accept )       !! one Monte Carlo step per spin
      call calc_E( spin, E )
      R_E( j ) = E
      R_spin(1:L,1:L,1:L, R_perm(j) ) = spin(1:L,1:L,1:L)   !! keep the state of the replica
!
      R_mti(j) = mti_mt64  ;  R_mt( : , j) = mt_mt64( : )   !! keep the mt19937_64 state (replica j)
    end do
!
    call replica_exchange( iter )      !! exchange R_perm(j) (or not) !! R_T(j) and R_E(j) corresponds to R_spin(:,:, R_perm(j) )
!
  END DO
!
! exchange MC, averaging observable data
  R_sum_observables( : , 1:R_num ) = 0.0d0                  !! reset summation E, E**2, M, M**2 etc
  R_accept_EXMC( 1:R_num ) = 0.0d0
  DO iter = 1, N_mcs_average                                !! MC steps for averaging
!$OMP parallel do default(none) &
!$OMP private( j, temperature, spin, E, M, accept ) &
!$OMP shared( L, R_num, R_T, R_spin, R_E, R_perm, R_mti, R_mt ) schedule(dynamic,1)
    do j = 1, R_num
      temperature = R_T(j)
      mti_mt64 = R_mti(j)  ;  mt_mt64( : ) = R_mt( : , j)   !! restore the mt19937_64 state (replica j)
!
      spin(1:L,1:L,1:L) = R_spin(1:L,1:L,1:L, R_perm(j) )   !! restore the state of the replica
      accept = 0.0d0
        CALL one_MC_step( temperature, spin, accept )       !! one Monte Carlo step per spin
      call calc_E( spin, E )
      M = DBLE( SUM( spin(1:L,1:L,1:L) ) )                  !!!!! call calc_M( spin, M )
      R_E( j ) = E
      CALL sum_observable_data( j, E, M, accept )
      R_spin(1:L,1:L,1:L, R_perm(j) ) = spin(1:L,1:L,1:L)   !! keep the state of the replica
!
      R_mti(j) = mti_mt64  ;  R_mt( : , j) = mt_mt64( : )   !! keep the mt19937_64 state (replica j)
    END DO
!
    call replica_exchange( iter )  !! exchange R_perm(j) (or not) !! R_T(j) and R_E(j) corresponds to R_spin(:,:,:, R_perm(j) )
!
  END DO
!
  CALL output_table( nf )
  close( nf )
!
  contains
!
!-----
SUBROUTINE one_MC_step( temperature, spin, accept )
  use,intrinsic :: iso_fortran_env
  use exchange_MC , only : L, N_spin
  use mt19937_64 , only : genrand64_real3
  IMPLICIT NONE
  real(REAL64), intent(in) :: temperature
  integer, intent (inout) :: spin(L,L,L)
  real(REAL64), intent(inout) :: accept
  integer :: i_spin, x,y,z
  real(REAL64) :: dE
!
  DO i_spin = 1, N_spin                                      !! one Monte Carlo step per spin
    x = int( L*genrand64_real3()+1 )  !! (x,y,z) of a trial spin
    y = int( L*genrand64_real3()+1 )
    z = int( L*genrand64_real3()+1 )
    call calc_dE( x,y,z, spin, dE )                          !! dE = E_{after trial} - E_{before trial}
    IF( genrand64_real3() <= exp( - dE/temperature ) ) THEN  !! accept ( Metropolis )
      spin(x,y,z) = - spin(x,y,z)
      accept = accept + 1.0d0
 !! ELSE                                                     !! reject
    END IF
  END DO
!
END SUBROUTINE one_MC_step
!
!-----
subroutine calc_dE( x,y,z, spin, dE )    !! calc dE = E_{after trial} - E_{before trial}
  use,intrinsic :: iso_fortran_env
  use exchange_MC , only : L, Jnn
  IMPLICIT NONE
  integer, intent(in) :: x,y,z, spin(L,L,L)
  real(REAL64), intent(out) :: dE
  integer :: xp1, xm1, yp1, ym1, zp1, zm1
      xm1 = x - 1  ;  IF(x == 1) xm1 = L
      xp1 = x + 1  ;  IF(x == L) xp1 = 1
      ym1 = y - 1  ;  IF(y == 1) ym1 = L
      yp1 = y + 1  ;  IF(y == L) yp1 = 1
      zm1 = z - 1  ;  IF(z == 1) zm1 = L
      zp1 = z + 1  ;  IF(z == L) zp1 = 1
  dE = Jnn * DBLE(2*spin(x,y,z)*( spin(xm1,y,z)+spin(xp1,y,z)+spin(x,ym1,z)+spin(x,yp1,z)+spin(x,y,zm1)+spin(x,y,zp1) ) )
end subroutine calc_dE
!
!-----
subroutine calc_E( spin, E )            !! calc  E = energy
  use,intrinsic :: iso_fortran_env
  use exchange_MC , only : L, Jnn
  IMPLICIT NONE
  integer, intent(in) :: spin(L,L,L)
  real(REAL64), intent(out) :: E
  integer :: x,y,z, xp1,yp1,zp1
  E = 0.0d0
  DO z = 1, L
    zp1 = z + 1  ;  IF(z == L) zp1 = 1
    DO y = 1, L
      yp1 = y + 1  ;  IF(y == L) yp1 = 1
      DO x = 1, L
        xp1 = x + 1  ;  IF(x == L) xp1 = 1
        E = E + DBLE( spin(x,y,z) * ( spin(x,yp1,z)+spin(xp1,y,z)+spin(x,y,zp1) ) )
      END DO
    END DO
  END DO
  E = - Jnn * E
end subroutine calc_E
!
!------
subroutine replica_exchange( iter )                  !! replica exchange
  use,intrinsic :: iso_fortran_env
  use exchange_MC , only : R_num, R_mti, R_mt, R_T, R_E, R_perm, R_accept_EXMC
  use mt19937_64 , only : genrand64_real3, mti_mt64 => mti, mt_mt64 => mt
  IMPLICIT NONE
  integer, intent(in) :: iter
  integer :: j_start, jj, j
  real(REAL64) :: Delta, xx
  mti_mt64 = R_mti( R_num + 1 )  ;  mt_mt64( : ) = R_mt( : , R_num + 1 )  !! restore the mt19937_64 state
!
  j_start = 1                                        !!  iter (MCstep) == 1, 3, 5, ...
  if( mod(iter,2) == 0 ) j_start = 2                 !!  iter (MCstep) == 2, 4, 6, ...
  do j = j_start, R_num - 1 , 2
    Delta = ( 1.0d0/R_T(j+1) - 1.0d0/R_T(j) )*( R_E(j) - R_E(j+1) )
    xx = exp( - Delta )
    if( xx < genrand64_real3() ) cycle               !! no exchange
               jj = R_perm( j   )                    !! exchange replicas j and (j+1)
    R_perm( j   ) = R_perm( j+1 )                    !! [R_T(j) and R_E(j)] corresponds to R_spin(:,:, R_perm(j) )
    R_perm( j+1 ) = jj
    R_accept_EXMC( j ) = R_accept_EXMC( j ) + 1.0d0  !! acceptance probability of exchange replicas at T = R_T(j)
  end do
!
  R_mti( R_num + 1 ) = mti_mt64  ;  R_mt( : , R_num + 1 ) = mt_mt64( : )  !! keep the mt19937_64 state
end subroutine replica_exchange
!
END PROGRAM main
!

