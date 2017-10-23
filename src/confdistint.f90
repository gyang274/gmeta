!subroutine CONFDISTINT(K,Thetahat,se,n,resample,B,len,m,band_pwr,shrink,smoothlist,eligid,distv,finalxi1,finalxi2,finalxi3)
subroutine CONFDISTINT(K,Thetahat,se,n,resample,B,len,m,band_pwr,shrink,smoothlist,eligids,distvs,finalxi1,finalxi2,finalxi3)
   implicit none
!------------------------------------   
! indicate variable type 
!-------------------------------------
   !------------------
   ! input
   !------------------
   integer:: resample, B, len, K,m
   double precision:: Thetahat(K), se(K), n(K), band_pwr

   
   !------------------
   ! output
   !------------------
   !double precision:: percentiles(3,2) 
   double precision:: shrink
   double precision:: smoothlist(len), p(len,len)
   
   !------------------
   ! intermediate 
   !------------------
   integer:: i,ip(K), j,bb,ii,jj,ISEED,eligcount,minid,maxid ,ipd(len**2),ml1,ml2,mr1,mr2,ipresamp(resample),ipB(B),eligids(len**2)
   integer:: seed
   double precision:: thord(K), y(K),xx
   double precision:: R,xi(K,resample),xim(resample), C, S,seq(len),ThetahatM, sstot,sswithin,THstar,Thetastar(K)
   double precision:: parray(len,len,B),Thetahatstar(K),xistar(K,resample),xistarm(resample)
   double precision:: xismootharr(len,len,resample),tempv
   double precision:: temp(len,len,B),psort(B,len,len),yp(B),tempM(len,len),tempMv(len**2)
   double precision:: dist(len,len),distv(len**2),distvs(len**2),ydist(len**2)
   double precision:: finalxi1(resample),finalxi2(resample),finalxi3(resample),lmin,rmin
   double precision:: thresh,lmax,rmax,lmean,rmean
   !double precision:: yxi(resample),yr(len**2) !for imsl
   double precision:: xsmarb(resample)
   
   double precision,allocatable:: eligm(:,:)
   integer,allocatable:: eligid(:),eligidsum(:),eligidm(:,:),ipelig(:)
   
   double precision,external::smooth, normal1
   integer,external:: idl,idr   
   

!------------------
!step1: preparation
!------------------
  !------------------
  !1.1 sort Thetahat & se
  !------------------
    ip=(/(i,i=1,K)/)
    !call sort_real(Thetahat, y, iperm=ip) !this is for imsl
    call shell_sort(Thetahat,K,ip) 
    !Thetahat = Thetahat(ip(1:K)) ! this is for imsl
    se = se(ip(1:K)) 
  !----------
  !1.2 get xi and xim
  !----------
    !ISEED = 123456   !for imsl
    !CALL RNSET (ISEED)   !for imsl
     do i=1,resample !xi is K*resample
       do j=1, K
         call nurnd(R) 
         R=R*se(j)+Thetahat(j)  
         !R=normal1(Thetahat(j),se(j))   !for earlier function
         xi(j,i)=R
       end do
       !call sort_real(xi(:,i), y)
       !xim(i)=y(m)
       ip=(/(ii,ii=1,K)/)
       y=xi(:,i)
       call shell_sort(y,K,ip)
       xim(i)=y(m)
     end do    
     
   !----------
   !1.3 C & S
   !----------
    C=sqrt(sum(n*(se**2))/dble(K))**(dble(1)-band_pwr)
    S=se(m)**band_pwr 
   !----------
   !1.4 smoothlist
   !---------- 
    seq=dble(log(0.01))
    do i=2,len
      seq(i)=seq(i-1)+(dble(log(3.0))-dble(log(0.01)))/dble(len-1)
    end do 
       
    smoothlist=exp(seq)*S*C

!------------------
!step2: tuning
!------------------
   !----------
   !2.1 shrink
   !----------  
    ThetahatM = sum(Thetahat)/dble(K)
    sstot=sum((Thetahat-ThetahatM)**2)/dble(K-1)
    sswithin=sum(se**2)/dble(K)
    shrink=min(dble(1.0),sswithin/sstot*(sqrt(sum(n*(se**2))/dble(K))/se(m))**0.1) 
    Thetastar=ThetahatM*shrink+Thetahat*(dble(1)-shrink)
    THstar=Thetastar(m)

   !----------
   !2.2 build prray :40 matrices of 10*10
   !----------
    parray=dble(0.0)
    
    do bb=1, B
      !--------
      !2.2.1 Thetahatstar
      !--------
      Thetahatstar=dble(0.0)
        do j=1,K
         !R=RNNOF()  !for imsl
         call nurnd(R) 
         R=Thetastar(j)+se(j)*R   
         !R=normal1(Thetastar(j),se(j))   
         Thetahatstar(j)=R
        end do
      !call sort_real(Thetahatstar, y)
      !Thetahatstar = y
      call shell_sort(Thetahatstar,K,ip) 
       
      !--------
      !2.2.2 xistar and xistarm
      !--------      
      do i=1,resample !xi is K*resample
        do j=1, K
          !R=RNNOF()   !for imsl
          call nurnd(R) 
          R=Thetahatstar(j)+se(j)*R   
          !R=normal1(Thetahatstar(j),se(j))   
          xistar(j,i)=R
        end do    
       !call sort_real(xistar(:,i), y)! these two lines are under imsl
       !xistarm(i)=y(m)   !this is xistar(i)!!!!!!!!!
       ip=(/(ii,ii=1,K)/)
       y=xistar(:,i)
       call shell_sort(y,K,ip)
       xistarm(i)=y(m)
      end do   
      !--------
      !2.2.3 xismootharr
      !--------
      xismootharr=dble(0.0) !200 matrices of 10*10
      do i=1,resample
        do ii=1,len !column
          do jj=1,len !row
            tempv=smooth(xistarm(i),xistar(:,i),smoothlist(jj),smoothlist(ii),K)
            xismootharr(jj,ii,i)=tempv
          end do
        end do
      end do
      !--------
      !2.2.4 compute p(,,bb)
      !-------- 
      do ii=1,len
        do jj=1,len
           !parray(jj,ii,bb)=sum((xismootharr(jj,ii,:)<=THstar)*(-1.0))/(1.0*resample)
           xsmarb=dble(0.0)
             do i=1,resample
               if (xismootharr(jj,ii,i)<=Thstar) then
               xsmarb(i)=dble(1.0)
               end if
             end do
           parray(jj,ii,bb)=sum(xsmarb)/dble(resample)
        end do
      end do   
      
    end do !end do bb, name it maybe
    
   !----------
   !2.3 now we find the smallest distance to uniform
   !----------   
    temp=0.0 !10*10*40
    
    do ii=1,len
      do jj =1,len
        !call sort_real(parray(jj,ii,:), yp)
        !psort(:,jj,ii)=yp !psort 40*10*10
        yp=parray(jj,ii,:)
        ipB=(/(i,i=1,B)/)
        call shell_sort(yp,B,ipB)
        psort(:,jj,ii)=yp 
      end do
    end do
    
    do bb=1,B
     temp(:,:,bb)=(psort(bb,:,:)-dble(bb/(B+1)))**2
    end do
    
    do ii=1,len
      do jj=1,len
        tempM(ii,jj)=sum(temp(ii,jj,:))/dble(B)
      end do
    end do
    
    tempMv=reshape(tempM,(/len**2/))
    distvs=tempMv-dble((/(i,i=1,len**2)/)/10000000)
    distv=distvs
    
    ipd=(/(i,i=1,len**2)/) 
    !call sort_real(distv, yr, iperm=ipd) !for imsl
    call shell_sort(distv,len**2,ipd)
    minid=ipd(1)
    lmin=smoothlist(idl(minid,len))
    rmin=smoothlist(idr(minid,len))

   !----------
   !2.4 next we use threshold to provide more possibility
   !----------  
    thresh=0.3596*(dble(B)**(-0.95))
    !eligcount=floor(sum((distv<thresh)*(-1.0)))
    eligids=(/(0,i=1,len**2)/) !record which is smaller than threshold
    eligcount=0
    do i=1,len**2
      if (distvs(i)<thresh) then
      eligcount=eligcount+1
      eligids(i)=i
      end if
    end do
    
    !write(*,*) "eligcount",eligcount

    if (eligcount>1)  then
      allocate(eligm(eligcount,2))
      allocate(eligid(eligcount)) 
      allocate(eligidm(eligcount,2)) 
      allocate(eligidsum(eligcount)) 
      allocate(ipelig(eligcount)) 
    
      ipd=(/(i,i=1,len**2)/)
      !call sort_real(distv, yr, iperm=ipd) !for imsl
      call shell_sort(distv,len**2,ipd)
      
      eligid=ipd(1:eligcount)
      
     
      do i=1,eligcount
        eligidm(i,1)=idl(eligid(i),len) !this is the index id, to find maxid
        eligidm(i,2)=idr(eligid(i),len)
        eligm(i,1)=smoothlist(idl(eligid(i),len))
        eligm(i,2)=smoothlist(idr(eligid(i),len))
      end do
      
      eligidsum=eligidm(:,1)+eligidm(:,2)
      
      !jj=ISMIN(eligcount,eligidsum*(-1.0),1)!for imsl      
      !maxid=eligid(jj)!for imsl
      ipelig=(/(i,i=1,eligcount)/)
      call shell_sort(dble(eligidsum),eligcount,ipelig) !convert the integer "eligidsum" to double precision
      jj=ipelig(eligcount) !find the max
      maxid=eligid(jj)
      
      lmax=smoothlist(idl(maxid,len))
      rmax=smoothlist(idr(maxid,len))
      
      lmean=sum(eligm(:,1))/eligcount
      rmean=sum(eligm(:,2))/eligcount
      
    else
    lmax=lmin
    lmean=lmin
    rmax=rmin
    rmean=rmin
    allocate(eligid(1))
    eligid=0 
    end if
    
!------------------
!step3: computing quantiles
!------------------
    do i=1,resample
    finalxi1(i)=smooth(xim(i),xi(:,i),lmin,rmin,K)
    finalxi2(i)=smooth(xim(i),xi(:,i),lmax,rmax,K)
    finalxi3(i)=smooth(xim(i),xi(:,i),lmean,rmean,K)
    end do

    
 end subroutine
        

!------------------------------------   
! smoothlist function
!------------------------------------
 double precision function smooth(cent, vec, l, r, K)
     integer K,j
     double precision vec(K), cent, l, r,S,vecc(K)
     vecc=dble(0.0)
     do j=1,K
      if ((vec(j)-cent)>=(-l) .and. (vec(j)-cent)<=r) then
      vecc(j)=dble(1.0)
      end if
     end do
     S=sum(vecc)
     smooth = dble(0)
     smooth = sum(vecc*vec)/S
     return
 end function smooth
   
!------------------
!locate id
!------------------
 integer function idl(x,len)
      integer x,len
      idl=mod(x,len)
       if (idl==0) then
       idl=10 
       end if
 end function idl
 
 integer function idr(x,len)
      integer x,len,idl
      idl=mod(x,len)
       if (idl==0) then
       idr=x/len
       else 
       idr=(x-idl)/len+1
       end if
 end function idr

!------------------
!shell sort
!------------------
subroutine shell_sort(a,ns,iperm)
  implicit none
  integer ns
  double precision :: a(ns) ! input
  integer is,js      ! do loop
  integer tempi      ! for exchange
  double precision temps
  integer ks         ! k 
  integer iperm(ns)
  ks=ns/2             ! initial
  do while( ks>0 )
    do is=ks+1,ns
      js=is-ks
      do while( js>0 )
      ! if a(j)>a(j+k) exchange
      ! a(j-k)\a(j) are new data
        if ( a(js) .gt. a(js+ks) ) then
          temps=a(js)
          a(js)=a(js+ks)
          a(js+ks)=temps
          tempi=iperm(js)
          iperm(js)=iperm(js+ks)
          iperm(js+ks)=tempi
          js=js-ks
        else
          exit ! a(j)<a(j+k) get out
        end if
   end do
 end do
    ks=ks/2 ! new k
 end do
  return
end subroutine
!------------------
!normal from Rmath
!------------------
subroutine nurnd(xx)
double precision normrnd, unifrnd, xx

call rndstart()
xx = normrnd()
call rndend()

return
end

!------------------
!returns random number between 0 - 1
!------------------
function ran1() 
implicit none
integer :: flag
double precision :: ran1
save flag
data flag /0/
if(flag==0) then
call random_seed(); flag = 1
endif
call random_number(ran1) ! built in fortran 90 random number function
return
end function ran1
!------------------
!returns a normal distribution
!------------------
function normal1(mean,sigma) 
implicit none
integer :: flag
double precision :: normal1, tmp, mean, sigma
double precision :: fac, gsave, rsq, r1, r2
double precision,external:: ran1
save flag, gsave
data flag /0/
if (flag.eq.0) then
rsq=2.0d0
do while(rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do
r1 = 2.0d0 * ran1() - 1.0d0
r2 = 2.0d0 * ran1() - 1.0d0
rsq = r1 * r1 + r2 * r2
enddo
fac = sqrt(-2.0d0*log(rsq)/rsq)
gsave = r1 * fac
tmp = r2 * fac
flag = 1
else
tmp = gsave
flag = 0
endif
normal1 = tmp * sigma + mean
return
end function normal1
