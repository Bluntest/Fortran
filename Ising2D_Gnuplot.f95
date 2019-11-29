MODULE CSEED
INTEGER :: ISEED
END MODULE CSEED
!
program Ising3D
USE CSEED
integer :: implicit none
INTEGER ::N=10000,Ni=20,Nj=20,i,j,l,n1,k,o
INTEGER :: time,STIME,seed
INTEGER, DIMENSION (9) :: TT
real *8 ::sum,E, M,M2,U,U2,C,Chi,Energy,Magnetization,kT,a(20,20),b,mu,x
real  :: RAND
b = 0. ! valor inicial do campo magnetico
x = -2. ! valor inicial do exchange J
do k=0, 10
  do o=1, 5
    if(x.ne.0) then
      Open(8, FILE='specific_heat.dat', Status='Unknown')
      Open(9, FILE='susceptibility.dat', Status='Unknown')
      Open(10, FILE='energy.dat', Status='Unknown')
      Open(11, FILE='magnetization.dat', Status='Unknown')
      !
      STIME = time()
      CALL gmtime(STIME,TT)
      ISEED = TT(6)+70*(TT(5)+12*(TT(4)+31*(TT(3)+23*(TT(2)+59*TT(1)))))
      IF (MOD(ISEED,2).EQ.0) ISEED = ISEED-1
      call srand (iseed)
      ! N=200 -Number of try
      ! N=2001
      ! Ni=20, Nj=20, Nk=20 - Lattice size
      ! a(20,20,20) - Spin configuration
      ! kT - Temperature
      write(*,*) '               kT                  C(Specific Heat)'
      write(*,*) '----------------------------------------------------'
      kT = 10.
      mu = 0.33
      1 continue
      !Running Averages
      ! Energy
	    U  = 0.0
      U2 = 0.0
      ! Magnetization
	    M  = 0.0
	    M2 = 0.0
      ! Set up the initial lattice configuration,
 	    do i = 1,Ni
        do j = 1,Nj
	        a(i,j) = 1.0
        end do
      end do
      ! Initial Energy and magnetization
	    Energy = 0.0
	    Magnetization = 0.0
	    do i = 1,Ni
	      do j = 1,Nj
	        sum=0.0 !sum of nearest neighbors
          ! Periodic Boundary Conditions
	        if(i>1 ) then
		        sum=sum+a(i-1,j)
	        endif
          if(i<Ni) then
		        sum=sum+a(i+1,j)
	        endif
	        if(j>1 ) then
            sum=sum+a(i,j-1)
	        endif
          if(j<Nj) then
		        sum=sum+a(i,j+1)
	        endif
	        Magnetization = Magnetization - a(i,j)
	        Energy=Energy - x*sum*a(i,j) - mu*b*a(i,j)
	      end do
	    end do
      !Correct the energy
	    Energy=Energy/2.0
	    do l = 1, N
        do n1 = 1,Ni*Nj
          !Pick the node to flip at random
          i = RAND()*Ni+1
	        j = RAND()*Nj+1
          ! Calculate	 change in energy	 when flipped
	         sum=0.0
	         if(i>1 )then
		         sum=sum+a(i-1,j)
	         endif
	         if(i<Ni) then
 		         sum=sum+a(i+1,j)
	         endif
	         if(j>1 ) then
		         sum=sum+a(i,j-1)
	         endif
	         if(j<Nj) then
	           sum=sum+a(i,j+1)
	         endif
	         E = 2.0*x*sum*a(i,j) + 2.0*mu*b*a(i,j) !energy change
           ! Metroplis criterium
	         if(( E<0.0 ).or.( exp(-E/kT)>RAND() )) then
             !If accept the move upgrade spin, total magnetization, and energy
		         a(i,j)=-a(i,j) !flip the spin
	           Energy = Energy + E
		         Magnetization = Magnetization - 2.0 * a(i,j)
	        endif
	      end do
        ! Increment averages
	      M = M + Magnetization
	      M2= M2+ Magnetization * Magnetization
	      U = U + Energy
	      U2= U2+ Energy * Energy
        ! Print out to files
	      if(modulo(l,N-1) == 0)then
	      U = U/ l/Ni/Nj !average energy
	      U2 = U2 / l/Ni/Nj/Ni/Nj
        C = ((U2-U*U))/kT/kT !specific heat
	      M = M / l /Ni/Nj !average magnetization
	      M2 = M2 / l /Ni/Nj/Ni/Nj
	      Chi = ((M2-M*M))/kT !susceptibility
        write(*,*) kT, C
	      write(8,*) kT, C ! print to file the specific heat
   	    write(9,*) kT, Chi !print to file the susceptibility
	      write(10,*) kT, U ! prin to file the average energy
	      write(11,*) kT, abs(M)! print to file the average magnetization
	      endif
      end do
      kT = kT+1.
      if(kT <= 100.0) go to 1
      close(8)
      close(9)
      close(10)
      close(11)
      42  format(a,f10.0,a,f10.0,a)
      open(1,file='png.plt',status='unknown')
        write(1,*) 'set terminal png'
        write(1,*) 'set grid'
        write(1,*) 'set xlabel "KT"'
        write(1,*) 'set ylabel "Energia"'
        !write(1,*) 'set yrange [:]'
        !write(1,*) 'set xrange [:]'
        write(1,42) 'set output "Energiab_',b,'J_',x,'.png"'
        write(1,42) 'plot "energy.dat" u 1:2 title "b_',b,'J_',x,'" w l lw 1.5 lt 1 lc 7'
        call system('gnuplot "png.plt"')
      close(1)
      open(1,file='png.plt',status='unknown')
        write(1,*) 'set terminal png'
        write(1,*) 'set grid'
        write(1,*) 'set xlabel "KT"'
        write(1,*) 'set ylabel "C (Calor especifico)"'
        !write(1,*) 'set yrange [:]'
        !write(1,*) 'set xrange [:]'
        write(1,42) 'set output "Calor-especifico-b_',b,'J_',x,'.png"'
        write(1,42) 'plot "specific_heat.dat" u 1:2 title "b_',b,'J_',x,'" w l lw 1.5 lt 1 lc 7'
        call system('gnuplot "png.plt"')
      close(1)
      open(1,file='png.plt',status='unknown')
        write(1,*) 'set terminal png'
        write(1,*) 'set grid'
        write(1,*) 'set xlabel "KT"'
        write(1,*) 'set ylabel "Susceptibilidade"'
        !write(1,*) 'set yrange [:]'
        !write(1,*) 'set xrange [:]'
        write(1,42) 'set output "Susceptibilidade-b_',b,'J_',x,'.png"'
        write(1,42) 'plot "susceptibility.dat" u 1:2 title "b_',b,'J_',x,'" w l lw 1.5 lt 1 lc 7'
        call system('gnuplot "png.plt"')
      close(1)
      open(1,file='png.plt',status='unknown')
        write(1,*) 'set terminal png'
        write(1,*) 'set grid'
        write(1,*) 'set xlabel "KT"'
        write(1,*) 'set ylabel "Magnetizacao"'
        !write(1,*) 'set yrange [:]'
        !write(1,*) 'set xrange [:]'
        write(1,42) 'set output "Magnetizacao-b_',b,'J_',x,'.png"'
        write(1,42) 'plot "magnetization.dat" u 1:2 title "b_',b,'J_',x,'" w l lw 1.5 lt 1 lc 7'
        call system('gnuplot "png.plt"')
      close(1)
    end if
    b = b + 10.
  end do
  x = x + 1.
end do
end program Ising3D
