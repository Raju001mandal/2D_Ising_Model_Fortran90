program ising_spin
implicit none


integer::i,j,k,d,ok,a,b,c,d2,n,step,mcs,ii,jj,p,q,z
real::num,stotal,m,energy,E,Eb,Ef,dE,u,h,s_av
integer::seed
integer,parameter::outunit=50
real,parameter::T=1.11,j_couple=1.0
integer,allocatable,dimension(:,:)::spin,test
character(len=30) ::fn

seed = 65452

print*,"give the dimension of the spin/test matrix(square)"
read*,d
print*,"give the number of monte carlo steps"
read*,step

allocate(spin(d,d),stat=ok)
!print*,"ok=",ok

i=p
j=q
energy=0.0
stotal=0.0
s_av=0.0
z=0
n=d*d
!================================initializing spin matrix===============================================

open(1,file="0.dat")

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
end do
write(1,fmt='(4g10.8)') (spin(j,i),j=1,d) 
end do
close (1)

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
    
   stotal=stotal+spin(i,j) 
   energy=energy-j_couple*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d2)))
  end do
end do   

m=stotal/real(n)
E=energy*0.5

print*,"Initial Magnetization= ",m
print*,"Initial Energy= ",E

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
 
 !=====================================================================================
  
  write(2,*)mcs,stotal/float(n),E/float(n)
! if(mcs.ge.250 )then 
  s_av=s_av+stotal/float(n)
  
! end if 
end do

s_av=s_av/step

s_av=abs(s_av)
print*,"average magnetization= ",s_av
close(2)
end program
        

