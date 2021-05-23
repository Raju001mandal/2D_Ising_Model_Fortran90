program ising_spin
implicit none


integer::i,j,k,d,ok,a,b,c,d2,n,step,mcs,ii,jj,p,q,z
real::num,stotal,m,energy,E,Eb,Ef,dE,u,h,s_av,m2,m2_av,e_av,dm2_av,e2,e2_av,sp_heat,m3,Enew
integer::seed
integer,parameter::outunit=50
real,parameter::T=1.11,j_couple=1.0
integer,allocatable,dimension(:,:)::spin,square
character(len=30) ::fn

seed = 65452

print*,"give the dimension of the spin/test matrix(square)"
read*,d
print*,"give the number of monte carlo steps"
read*,step

allocate(spin(d,d),stat=ok)
allocate(square(d,d),stat=ok)
!print*,"ok=",ok

i=p
j=q
energy=0.0
stotal=0.0
m2=0.0
s_av=0.0
m2_av=0.0
e_av=0.0
e2=0.0
e2_av=0.0
dm2_av=0.0

z=0
n=d*d  !totatl lattice sites
m2=0.0 !m2 will calculate average of spin square
!================================initializing spin matrix===============================================

open(1,file="0.dat")
open(11,file="sq.dat")

do i=1,d
  do j=1,d
    call random_number(num)
     
!  spin(j,i)=+1
!     spin(i,j)=-1
     
  if(num<0.5)then
   spin(j,i)=-1

   else
   spin(j,i)=+1

   end if
square(j,i)=spin(j,i)*spin(j,i)
end do
write(1,fmt='(4g10.8)') (spin(j,i),j=1,d) 
write(11,fmt='(4g10.8)') (square(j,i),j=1,d) 
end do
close (1)
close(11)
!================================calculation of initial magnetization and energy============================

do i=1,d
  do j=1,d
  
!!========neighbours============================= 
    a=i+1
    b=i-1
    c=j+1
    d2=j-1
!!========periodic boundary condition=============
   if(i==1)then
   b=d
   end if
   
   if(i==d)then
   a=1
   end if
   
   if(j==1)then
   d2=d
   end if
   
   if(j==d)then
   c=1
   end if
    
   stotal=stotal+spin(i,j) !summing over all the spins
   m2=m2+square(i,j)       !summing over square of spins
   energy=energy-j_couple*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d2)))
  end do
end do   

m=stotal/real(n)
m2=m2/real(n)      ! averaging of spin squares
E=energy*0.5
e2=E**2
print*,"Initial Magnetization= ",m
print*,"Initial sum of spin square m2 = ",m2
print*,"Initial EnergyE,E2= ",E,e2

!===========================monte carlo step to reach equilibrium==================================

open(2,file="ising.dat")

do mcs=1,step

  do ii=1,d
     do jj=1,d
     !=========choosing a random lattice site====================
      call random_number(num)
      i=int(num*float(d))+1
      call random_number(num)
      j=int(num*float(d))+1
      
      !!========neighbours============================= 
       a=i+1
       b=i-1
       c=j+1
       d2=j-1
     !=============periodic boundary condition ===================
      if(i==1)then
      b=d
      end if
   
      if(i==d)then
      a=1
      end if
   
      if(j==1)then
      d2=d
      end if
   
      if(j==d)then
      c=1
      end if  
     !==========================energy calculation before spin flip=========================================
      Eb=-j_couple*float(spin(i,j)*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d2)))
     !==========================spin flip=========================================
      spin(i,j)=-spin(i,j)
 !     square(i,j)=spin(i,j)*spin(i,j)
    
     !==========================energy calculation after spin flip================
     Ef=-j_couple*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d2)))    
     dE=Ef-Eb
     !=========accepting the flip and updating the energy and magnetization=================
      if(dE<=0.0)then
        E=E+dE
      
        stotal=stotal+2*float(spin(i,j)) !change in spin (1-(-1))
  
      else
        u=exp(-dE/(T))
        call random_number(h)
        if(h<u)then
         E=E+dE
        
         stotal=stotal+2*float(spin(i,j))
 
        else
         spin(i,j)=-spin(i,j)
 
        end if
      end if
      
     end do
   end do
 !============================for visualization of microstates=========================================================
 
if(mod(mcs,10)==0)then
  write(fn,fmt='(i0,a)') mcs,'.dat'
  open(unit=outunit,file=fn,form='formatted')
    do i=1,d
     write(outunit,*)(spin(i,j),j=1,d) 
    end do
    close(outunit) 
 end if  
!===============================================================================================
  Enew=E/float(n)
  e2=Enew**2
   
  m3=stotal/float(n)
  m2=(m3)**2
 !=====================================================================================
  
write(2,*)mcs,m3,E,E/float(n),e2,m2

  s_av=s_av+m3
  m2_av=m2_av+m2
  e_av=e_av+Enew
  e2_av=e2_av+e2

end do

print*,"m2= ",m2_av

s_av=s_av/real(step) !calculation of average m (ensemble average for per spin magnetization.)
s_av=abs(s_av)
m2_av=m2_av/real(step)
dm2_av=abs((m2_av-(s_av)**2))
e_av=e_av/real(step)
e2_av=e2_av/real(step)
sp_heat=abs((e2_av-(e_av)**2)/(T**2))
print*,"average magnetization= ",s_av
print*,"m2= ",m2_av
print*,"average energy= ",e_av
print*,"average of square energy= ",e2_av
print*,"specific heat= ",sp_heat
print*,"Fluctuation= ",dm2_av
close(2)
end program
        

