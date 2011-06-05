
!==========================================================
!
!      Latent class mixed model for curvilinear outcome
!           using  Beta CDF or quadratic I-Splines
!
!        Cecile Proust, Helene Jacqmin-Gadda
!
!
!       Corresponding author :
!       Cecile Proust, INSERM U897, ISPED,
!       146  rue L\'eo Saignat,
!       33076 Bordeaux cedex, France.
!       Tel: (33) 5 57 57 45 79; Fax: (33) 5 56 24 00 81;
!       e-mail : cecile.proust@isped.u-bordeaux2.fr
!
!                                       21/12/2010
!===========================================================
! - Version fortran 90
!






!----------------------------------------------------------
!
!- Module COMMUN avec les donnees dynamiques
!
!----------------------------------------------------------




      module communc

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg  &
      ,nprob,nvarprob,maxmes,nobs,ntrtot,nrisq,nvarxevt,nef &
      ,idlink
      double precision,dimension(:),allocatable,save::Y
      double precision,dimension(:,:),allocatable,save ::X
      integer,dimension(:),allocatable,save ::idea,idg,idprob
      integer,dimension(:),allocatable,save :: nmes,prior
      double precision,save :: minY,maxY,epsY
      double precision,dimension(:),allocatable,save::zitr
      double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2

      end module communc


!----------------------------------------------------------
!
!     INTERFACE TYPEc
!
!----------------------------------------------------------

      module typec

      interface verif1c
      subroutine marq98c(b,m,ni,v,rl,ier,istop,ca,cb,dd)
         integer,intent(in) :: m
         integer,intent(inout)::ni,ier,istop
         double precision,dimension(m*(m+3)/2),intent(out)::v
         double precision,intent(out)::rl
         double precision,dimension(m),intent(inout)::b
	 double precision,intent(inout)::ca,cb,dd
      end subroutine marq98c

      subroutine derivac(b,m,v,rl)
        integer,intent(in)::m
        double precision,intent(inout)::rl
        double precision,dimension(m),intent(in)::b
        double precision,dimension((m*(m+3)/2)),intent(out)::v
      end subroutine derivac

      subroutine searpasc(vw,step,b,bh,m,delta,fim)
        integer,intent(in)::m
        double precision,dimension(m),intent(in)::b
        double precision,dimension(m),intent(inout)::bh,delta
        double precision,intent(inout)::vw,fim,step
      end subroutine searpasc

      subroutine dmfsdc(a,n,eps,ier)
        integer,intent(in)::n
        integer,intent(inout)::ier
        double precision,intent(inout)::eps
        double precision,dimension(n*(n+1)/2),intent(inout)::A
      end subroutine dmfsdc

      subroutine valfpac(vw,fi,b,bk,m,delta)
        integer,intent(in)::m
	double precision,intent(in)::vw
        double precision,dimension(m),intent(in)::b,delta
        double precision,dimension(m),intent(out)::bk
        double precision,intent(out)::fi
      end subroutine valfpac

      subroutine dmaxtc(maxt,delta,m)
        integer,intent(in)::m
        double precision,dimension(m),intent(in)::delta
        double precision,intent(out)::maxt
      end subroutine dmaxtc
      end interface verif1c

      interface verif2c
      subroutine dsinvc(A,N,EPS,IER,DET)
        integer,intent(in)::n
        integer,intent(inout)::ier
        double precision,intent(inout)::eps
        double precision,intent(inout),optional::det
        double precision,dimension(n*(n+1)/2),intent(inout)::A
      end subroutine dsinvc

      subroutine dcholec(a,k,nq,idpos)
      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a
      end subroutine dcholec
      end interface verif2c

      end module typec



!----------------------------------------------------------
!
!     MODULE PARAMETERS
!
! Derniere mise a jour : 09/02/2011
!-----------------------------------------------------------

      module parametersc
          double precision,save::epsa,epsb,epsd
          integer,save::maxiter
      end module parametersc

!-------------------------------------------------------------
!
!          MODULE OPTIM avec MARQ98
!
!-------------------------------------------------------------



      module optimc
      implicit none
! -Interface permettant la verification des type des arguments
      interface verif1c
        module procedure marq98c,derivac,searpasc,dmfsdc,valfpac
      end interface verif1c

      interface verif2c
        module procedure dsinvc,dcholec,dmaxtc
      end interface verif2c

      CONTAINS
!-------------------------------------------------------------
!                   MARQ98
!-------------------------------------------------------------


      subroutine marq98c(b,m,ni,v,rl,ier,istop,ca,cb,dd)

!
!  fu = matrice des derivees secondes et premieres
!
!  istop: raison de l'arret
!  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
!  2: nb max d'iterations atteints
!  4: Erreur

      use parametersc

      IMPLICIT NONE
!   variables globales
      integer,intent(in) :: m
      integer,intent(inout)::ni,ier,istop
      double precision,dimension(m*(m+3)/2),intent(out)::v
      double precision,intent(out)::rl
      double precision,dimension(m),intent(inout)::b
      double precision,intent(inout)::ca,cb,dd

!   variables locales
      integer::nql,ii,nfmax,idpos,ncount,id,jd,m1,j,i,ij
      double precision,dimension(m*(m+3)/2)::fu
      double precision,dimension(m)::delta,b1,bh
      double precision::da,dm,ga,tr
      double precision::GHG,funcpac,det,step,eps,vw,fi,maxt, &
      z,rl1,th,ep


      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-7!1.d-6
      nfmax=m*(m+1)/2
      ca=epsa+1.d0
      cb=epsb+1.d0
      rl1=-1.d+10
      ni=0
      istop=0
      da=0.01d0
      dm=5.d0
      nql=1
      m1=m*(m+1)/2
      ep=1.d-20

      Main:Do
!	write(*,*)'iteration',ni
	call derivac(b,m,v,rl)
!	write(*,*)'vrais',rl
        rl1=rl
        dd = 0.d0
        fu=0.D0
        do i=1,m
           do j=i,m
              ij=(j-1)*j/2+i
              fu(ij)=v(ij)
           end do
        end do

        call dsinvc(fu,m,ep,ier,det)
        if (ier.eq.-1) then
           dd=epsd+1.d0
        else
           GHG = 0.d0
           do i=1,m
              do j=1,m
                 if(j.ge.i) then
                    ij=(j-1)*j/2+i
                 else
                    ij=(i-1)*i/2+j
                 end if
                 GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
	      end do
           end do
           dd=GHG/dble(m)
        end if

        if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) exit main
        tr=0.d0
        do i=1,m
           ii=i*(i+1)/2
           tr=tr+dabs(v(ii))
        end do
        tr=tr/dble(m)

        ncount=0
        ga=0.01d0
 400    do i=1,nfmax+m
           fu(i)=v(i)
        end do
        do i=1,m
           ii=i*(i+1)/2
           if (v(ii).ne.0) then
              fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
           else
              fu(ii)=da*ga*tr
           endif
        end do
        call dcholec(fu,m,nql,idpos)
        if (idpos.ne.0) then
           ncount=ncount+1
           if (ncount.le.3.or.ga.ge.1.d0) then
              da=da*dm
           else
              ga=ga*dm
              if (ga.gt.1.d0) ga=1.d0
           endif
           goto 400
        else
            do i=1,m
               delta(i)=fu(nfmax+i)
               b1(i)=b(i)+delta(i)
            end do
            rl=funcpac(b1,m,id,z,jd,z)
            if (rl1.lt.rl) then
               if(da.lt.eps) then
                  da=eps
               else
                  da=da/(dm+2.d0)
               endif
               goto 800
            endif
         endif
!      write(6,*) 'loglikelihood not improved '
         call dmaxtc(maxt,delta,m)
         if(maxt.eq.0.D0) then
            vw=th
         else
            call dmaxtc(maxt,delta,m)
            vw=th/maxt
         endif
         step=dlog(1.5d0)
!      write(*,*) 'searpas'
         call searpasc(vw,step,b,bh,m,delta,fi)
         rl=-fi
         if(rl.eq.-1.D9) then
               istop=4
               goto 110
          end if

         do i=1,m
            delta(i)=vw*delta(i)
         end do
         da=(dm-3.d0)*da

 800     cb=dabs(rl1-rl)
         ca=0.d0
         do i=1,m
            ca=ca+delta(i)*delta(i)
         end do
!         write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
         do i=1,m
            b(i)=b(i)+delta(i)
         end do

         ni=ni+1
         if (ni.ge.maxiter) then
            istop=2
!            write(6,*) 'maximum number of iteration reached'
            goto 110
         end if
      End do Main
      v=0.D0
      v(1:m*(m+1)/2)=fu(1:m*(m+1)/2)
      istop=1

 110   continue
       return
       end subroutine marq98c

!------------------------------------------------------------
!                          DERIVA
!------------------------------------------------------------

      subroutine derivac(b,m,v,rl)

      implicit none

      integer,intent(in)::m
      double precision,intent(inout)::rl
      double precision,dimension(m),intent(in)::b
      double precision,dimension((m*(m+3)/2)),intent(out)::v
      double precision,dimension(m)::fcith
      integer ::i0,m1,ll,i,k,j
      double precision::funcpac,thn,th,z,vl,temp,thi,thj
!
!     v:matrice d'information+score
!     calcul de la derivee premiere par "central difference"
!     calcul des derivees secondes par "forward difference"
!
      z=0.d0
      i0=0

      rl=funcpac(b,m,i0,z,i0,z)

      if(rl.eq.-1.d9) then
         goto 123
      end if

      do i=1,m
         th=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
         fcith(i)=funcpac(b,m,i,th,i0,z)
         if(fcith(i).eq.-1.d9) then
            rl=-1.d9
            goto 123
         end if
      end do

      k=0
      m1=m*(m+1)/2
      ll=m1
      Main:do i=1,m
         ll=ll+1
         thn=-DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
         temp=funcpac( b,m, i,thn,i0,z)
         if(temp.eq.-1.d9) then
            rl=-1.d9
            exit Main
         end if
         vl=(fcith(i)-temp)/(2.d0*(-thn))
         v(ll)=vl
         do j=1,i
            k=k+1

            thi=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
            thj=DMAX1(1.d-7, 1.d-4 * DABS(b(j)))

            temp=funcpac(b,m,i,thi,j,thj)
            if(temp.eq.-1.d9) then
               rl=-1.d9
               exit Main
            end if
            v(k)=-(temp-fcith(j)-fcith(i)+rl)/(thi*thj)
         end do
      end do Main
 123   continue

      return
      end subroutine derivac

!------------------------------------------------------------
!                        SEARPAS
!------------------------------------------------------------


      subroutine searpasc(vw,step,b,bh,m,delta,fim)
!
!  MINIMISATION UNIDIMENSIONNELLE
!
      implicit none

      integer,intent(in)::m
      double precision,dimension(m),intent(in)::b
      double precision,intent(inout)::vw
      double precision,dimension(m),intent(inout)::bh,delta
      double precision,intent(inout)::fim,step
      double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3
      integer::i

       vlw1=dlog(vw)
       vlw2=vlw1+step
       call valfpac(vlw1,fi1,b,bh,m,delta)
       call valfpac(vlw2,fi2,b,bh,m,delta)

       if(fi2.ge.fi1) then
	  vlw3=vlw2
	  vlw2=vlw1
	  fi3=fi2
	  fi2=fi1
	  step=-step

          vlw1=vlw2+step
          call valfpac(vlw1,fi1,b,bh,m,delta)
          if(fi1.gt.fi2) goto 50
       else
          vlw=vlw1
          vlw1=vlw2
          vlw2=vlw
          fim=fi1
          fi1=fi2
          fi2=fim
       end if

       do i=1,40
          vlw3=vlw2
          vlw2=vlw1
          fi3=fi2
          fi2=fi1

          vlw1=vlw2+step
          call valfpac(vlw1,fi1,b,bh,m,delta)
          if(fi1.gt.fi2) goto 50
          if(fi1.eq.fi2) then
             fim=fi2
             vm=vlw2
             goto 100
          end if
       end do
!
!  PHASE 2 APPROXIMATION PAR QUADRIQUE
!
50     continue
!
!  CALCUL MINIMUM QUADRIQUE
!
      vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))
      call valfpac(vm,fim,b,bh,m,delta)
      if(fim.le.fi2) goto 100
      vm=vlw2
      fim=fi2
100   continue
      vw=dexp(vm)

      return

      end subroutine searpasc

!------------------------------------------------------------
!                         DCHOLE
!------------------------------------------------------------

      subroutine dcholec(a,k,nq,idpos)

      implicit none

      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a

      integer::i,ii,i1,i2,i3,m,j,k2,jmk
      integer::ijm,irm,jji,jjj,l,jj,iil,jjl,il
      integer,dimension(k)::is
      double precision ::term,xn,diag,p
      equivalence (term,xn)


!      ss programme de resolution d'un systeme lineaire symetrique
!
!       k ordre du systeme /
!       nq nombre de seconds membres
!
!       en sortie les seconds membres sont remplaces par les solutions
!       correspondantes
!

      i2=0
      ii=0
      idpos=0
      k2=k+nq
!     calcul des elements de la matrice
      do i=1,k
         ii=i*(i+1)/2
!       elements diagonaux
         diag=a(ii)
         i1=ii-i
         if(i-1.ne.0) goto 1
         if(i-1.eq.0) goto 4
1        i2=i-1
         do l=1,i2
             m=i1+l
             p=a(m)
             p=p*p
             if(is(l).lt.0) goto 2
             if(is(l).ge.0) goto 3
2            p=-p
3            diag=diag-p
         end do

4        if(diag.lt.0) goto 5
         if(diag.eq.0) goto 50
         if(diag.gt.0) goto 6
5        is(i)=-1
         idpos=idpos+1
         diag=-dsqrt(-diag)
         a(ii)=-diag
         goto 7
6        is(i)=1
         diag=dsqrt(diag)
         a(ii)=diag
!       elements non diagonaux
7        i3=i+1
         do j=i3,k2
            jj=j*(j-1)/2+i
            jmk=j-k-1
            if(jmk.le.0) goto 9
            if(jmk.gt.0) goto 8
8           jj=jj-jmk*(jmk+1)/2
9           term=a(jj)
            if(i-1.ne.0) goto 10
            if(i-1.eq.0) goto 13
10          do l=1,i2
               iil=ii-l
               jjl=jj-l
               p=a(iil)*a(jjl)
               il=i-l
               if(is(il).lt.0) goto 11
               if(is(il).ge.0) goto 12
11             p=-p
12             term=term-p
            end do
13            a(jj)=term/diag
	   end do
      end do

!       calcul des solutions
      jj=ii-k+1
      do l=1,nq
         jj=jj+k
         i=k-1
14       jji=jj+i
         xn=a(jji)
         if(i-k+1.lt.0) goto 20
         if(i-k+1.ge.0) goto 22
20       j=k-1
21       jjj=jj+j
         ijm=i+1+j*(j+1)/2
         xn=xn-a(jjj)*a(ijm)
         if(j-i-1.le.0) goto 22
         if(j-i-1.gt.0) goto 30
30       j=j-1
         goto 21
22       irm=(i+1)*(i+2)/2
         a(jji)=xn/a(irm)
         if(i.le.0) cycle
         if(i.gt.0) goto 40
40       i=i-1
         go to 14
      end do
50    continue
      return
      end subroutine dcholec


      subroutine dmfsdc(a,n,eps,ier)
!
!   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
!   MATRICE = TRANSPOSEE(T)*T
!   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
!            PAR COLONNE DE LA METRICE A FACTORISER
!   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
!
!   SUBROUTINE APPELE PAR DSINV
!
!   N : DIM. MATRICE
!   EPS : SEUIL DE TOLERANCE
!   IER = 0 PAS D'ERREUR
!   IER = -1 ERREUR
!   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(inout)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

!
!   TEST ON WRONG INPUT PARAMETER N
!
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
!
!   INITIALIZE DIAGONAL-LOOP
!
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
!
!   CALCULATE TOLERANCE
!
          tol=dabs(eps*sngl(A(kpiv)))
!
!   START FACTORIZATION-LOOP OVER K-TH ROW
!
         do i=k,n
	    dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
!
!   START INNER LOOP
!
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
	       dsum=dsum+A(lanf)*A(lind)
            end do

!
!   END OF INNEF LOOP
!
!   TRANSFORM ELEMENT A(IND)
!
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
!
!   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
!
5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
!
!   COMPUTE PIVOT ELEMENT
!
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
!
!   CALCULATE TERMS IN ROW
!
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
	 ind=ind+i
      end do

!
!   END OF DIAGONAL-LOOP
!
      return
12    ier=-1
      return

      end subroutine dmfsdc


!------------------------------------------------------------
!                            DSINV
!------------------------------------------------------------


      subroutine dsinvc(A,N,EPS,IER,DET)

!
!     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
!
!     MATRICE = TRANSPOSEE(T)*T
!     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
!
!     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
!         STOCKEE COLONNE PAR COLONNE
!     DIM. MATRICE A INVERSER = N
!     DIM. TABLEAU A = N*(N+1)/2
!
!     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
!           COMME NUL
!
!     IER : CODE D'ERREUR
!         IER=0 PAS D'ERREUR
!         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
!         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(inout)::eps
      double precision,intent(inout),optional::det
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision::din,work
      integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf

!
!     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
!     A=TRANSPOSE(T) * T
!

      call dmfsdc(A,n,eps,ier)



      if (ier.lt.0) goto 9
      if (ier.ge.0) det=0.d0
!
!     INVERT UPPER TRIANGULAR MATRIX T
!     PREPARE INVERSION-LOOP
!
!
! calcul du log du determinant

      do i=1,n
         det=det+dlog(A(i*(i+1)/2))
      end do
      det=2*det
      ipiv=n*(n+1)/2
      ind=ipiv
!
!     INITIALIZE INVERSION-LOOP
!
      do i=1,n
         din=1.d0/A(ipiv)
         A(ipiv)=din
         min=n
         kend=i-1
         lanf=n-kend
         if (kend.le.0) goto 5
         if (kend.gt.0) j=ind
!
!     INITIALIZE ROW-LOOP
!
         do k=1,kend
	    work=0.d0
	    min=min-1
	    lhor=ipiv
	    lver=j
!
!     START INNER LOOP
!
            do l=lanf,min
	        lver=lver+1
		lhor=lhor+l
                work=work+A(lver)*A(lhor)
	    end do
!
!     END OF INNER LOOP
!
            A(j)=-work*din
            j=j-min
	 end do

!
!     END OF ROW-LOOP
!
5        ipiv=ipiv-min
         ind=ind-1
      end do

!
!     END OF INVERSION-LOOP
!
!     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
!     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
!     INITIALIZE MULTIPLICATION-LOOP
!
      do i=1,n
         ipiv=ipiv+i
	 j=ipiv
!
!     INITIALIZE ROW-LOOP
!
	 do k=i,n
	    work=0.d0
	    lhor=j
!
!     START INNER LOOP
!
            do l=k,n
	        lver=lhor+k-i
		work=work+A(lhor)*A(lver)
   		lhor=lhor+l
            end do
!
!     END OF INNER LOOP
!
            A(j)=work
            j=j+k
	 end do
      end do

!
!     END OF ROW-AND MULTIPLICATION-LOOP
!
9     return
      end subroutine dsinvc

!------------------------------------------------------------
!                          VALFPA
!------------------------------------------------------------

        subroutine valfpac(vw,fi,b,bk,m,delta)

        implicit none

        integer,intent(in)::m
        double precision,dimension(m),intent(in)::b,delta
        double precision,dimension(m),intent(out)::bk
        double precision,intent(out)::fi
	double precision::vw,funcpac,z
	integer::i0,i

         z=0.d0
         i0=1
         do i=1,m
            bk(i)=b(i)+dexp(vw)*delta(i)
	 end do
         fi=-funcpac(bk,m,i0,z,i0,z)

         return

         end subroutine valfpac

!------------------------------------------------------------
!                            MAXT
!------------------------------------------------------------


      subroutine dmaxtc(maxt,delta,m)

      implicit none

       integer,intent(in)::m
       double precision,dimension(m),intent(in)::delta
       double precision,intent(out)::maxt
       integer::i

       maxt=Dabs(delta(1))
       do i=2,m
         if(Dabs(delta(i)).gt.maxt)then
	    maxt=Dabs(delta(i))
	 end if
       end do

       return
       end subroutine dmaxtc

      end module optimc






!===========================================================
!
      !SUBROUTINES
!
!===========================================================


!-----------------------------------------------------------
!                    SUBROUTINE PRINCIPALE: HETMIXCONT
!------------------------------------------------------------



      subroutine hetmixCont(Y0,X0,Prior0,idprob0,idea0,idg0,ns0,ng0,nv0,nobs0 &
          ,nea0,nmes0,idiag0,nwg0,npm0,b,Vopt,vrais,ni,istop,gconv,ppi0,resid_m &
          ,resid_ss,pred_m_g,pred_ss_g,pred_RE,convB,convL,convG,maxiter0 &
          ,epsY0,idlink0,nbzitr0,zitr0,marker,transfY,nsim0)

      use parametersc
      use communc
      use optimc

      IMPLICIT NONE

	!D�claration des variables en entr�e
      integer,intent(in)::nv0,maxiter0
      double precision,intent(in)::epsY0
      integer,intent(in)::idlink0,nbzitr0
      double precision,dimension(nbzitr0),intent(in)::zitr0

      integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,npm0,nea0,nsim0
      integer, dimension(nv0),intent(in)::idea0,idg0,idprob0
      integer, dimension(ns0),intent(in)::nmes0,prior0
      double precision,dimension(nobs0),intent(in)::Y0
      double precision,dimension(nobs0*nv0),intent(in)::X0
      double precision,intent(in)::convB,convL,convG
	!D�claration des variable en entr�e et sortie
      double precision, dimension(npm0), intent(inout) :: b
	!D�claration des variables en sortie
      double precision,intent(out)::vrais
      double precision,dimension(3),intent(out)::gconv
      double precision,dimension(ns0*ng0),intent(out)::ppi0
      double precision,dimension(nobs0),intent(out)::resid_m,resid_ss
      double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
      double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
      double precision,dimension(ns0*nea0),intent(out)::pred_RE
      double precision,dimension(nsim0),intent(out)::marker,transfY
      double precision,dimension(npm0*(npm0+1)/2),intent(out)::Vopt
      integer, intent(out)::ni,istop
	!Variables locales
      integer::jtemp,i,g,j,ij,npm,ier,k,ktemp,ig,nmestot,it
      double precision::eps,ca,cb,dd
      double precision,dimension(ns0,ng0)::PPI
      double precision,dimension(npm0)::mvc
      double precision,dimension(npm0*(npm0+3)/2)::V


! sorties initialisees

       ppi0=0.d0
       Vopt=0.d0
       gconv=0.d0
       pred_ss_g=0.d0
       pred_m_g=0.d0
       pred_RE=0.d0
       marker=0.d0
       transfY=0.d0
       resid_m=0.d0
       resid_ss=0.d0
       vrais=0.d0
       ni=0





! en prevision de l'extension au conjoint
      nrisq=0
      nvarxevt=0
! fin en prevision

      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do

      epsa=convB
      epsb=convL
      epsd=convG
      maxiter=maxiter0

      minY=zitr0(1)
      maxY=zitr0(nbzitr0)
      epsY=epsY0
      idlink=idlink0
      if (idlink.eq.0) ntrtot=2
      if (idlink.eq.1) ntrtot=4
      if (idlink.eq.3) ntrtot=zitr0(nbzitr0)-zitr0(1)
      if (idlink.eq.2) then
         ntrtot=nbzitr0+2
         allocate(zitr(-1:(ntrtot)))

         allocate(mm(nobs0),mm1(nobs0),mm2(nobs0),im(nobs0),im1(nobs0),im2(nobs0))

         zitr(1:nbzitr0)=zitr0(1:nbzitr0)
         zitr(-1)=zitr(1)
         zitr(0)=zitr(1)
         zitr(ntrtot-1)=zitr(ntrtot-2)
         zitr(ntrtot)=zitr(ntrtot-1)
	else
         allocate(zitr(1))
         allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
      end if


      allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0) &
      ,idea(nv0),idg(nv0),nmes(ns0),prior(ns0))

      eps=1.d-20

! enrigstrement pour les modules
      ns=ns0
      ng=ng0
      nv=nv0
      nobs=nobs0
      if (nwg0.eq.0) then
         nwg=0
      else
         nwg=ng-1
      end if

      idiag=idiag0

      nmes=0
      Y=0.d0
      X=0.d0
      idprob=0
      idea=0
      idg=0
      nmestot=0
      ktemp=0
      do k=1,nv
         idprob(k)=idprob0(k)
         idea(k)=idea0(k)
         idg(k)=idg0(k)

         jtemp=0
         it=0
         DO i=1,ns
            if (k.eq.1) then
               nmes(i)=nmes0(i)
               prior(i)=prior0(i)
               do j=1,nmes(i)
                  nmestot=nmestot+1
                  jtemp=jtemp+1
                  Y(jtemp)=Y0(jtemp)
               end do
            end if

            do j=1,nmes(i)
               ktemp=ktemp+1
		       it=it+1
               X(it,k)=X0(ktemp)
            end do
         end do
      end do

! creation des parametres

      nea=0
      ncg=0
      ncssg=0
      nprob=ng-1
      nvarprob=min(ng-1,1)
      do k=1,nv
         if (idg(k).eq.1) then
            ncssg=ncssg+1      ! nb var. sans melange
         else if (idg(k).eq.2) then
            ncg=ncg+1      ! nb var. sans melange
         end if
         nea=nea+idea(k)
         nprob=nprob+(idprob(k))*(ng-1)
         nvarprob=nvarprob+idprob(k)
      end do

      if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
         istop=12
         go to 1589
      end if


!  nb effets fixes = nb effets fixes sans melange
!                  + ng fois le nb de var dans melange


      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

      nef=nprob+ncssg+ncg*ng+nvarxevt+nrisq-1
      npm=nef+nvc+nwg+ntrtot



      if (idiag.eq.1) then
         DO j=1,nvc
            B(nef+j)=dsqrt(abs(B(nef+j)))
         END DO
      end if

! si idiag=0, on met dans le vecteur des parms, les parms
! de la transformee de Cholesky

      if (idiag.eq.0) then

         DO j=1,nvc
            mvc(j)=B(nef+j)
         END DO

         CALL dmfsdc(mvc,nea,EPS,IER)
         DO j=1,nvc
            B(nef+j)=mvc(j)
         END DO
      end if
      if (nwg.gt.0) then
         do i=1,nwg
            B(nef+nvc+i)=abs(B(nef+nvc+i))
         end do
      end if


      if (idlink.eq.2) then
         call design_splines(ier)
         if (ier.eq.-1) then
            istop=9
            go to 1589
         end if
      end if

! lancement de l'optimisation

      IF (npm.eq.1) then
         istop=10
         go to 1589
      else
         ca=0.d0
         cb=0.d0
         dd=0.d0

         call marq98c(b,npm,ni,V,vrais,ier,istop,ca,cb,dd)

!         write(*,*)
!         write(*,*)'    FIN OPTIMISATION  ..... '
!         write(*,*)'istop',istop,'vrais',vrais



         gconv=0.d0
         gconv(1)=ca
         gconv(2)=cb
         gconv(3)=dd
         vopt(1:(npm*(npm+1)/2))=V(1:(npm*(npm+1)/2))

         do k=1,nwg
            b(nef+nvc+k)=abs(b(nef+nvc+k))
         end do

! probas posteriori

!      write(*,*)'avant postprobc'


         if (istop.eq.1) then
            if (ng.gt.1) then
               call postprobc(B,npm,PPI)
            end if



!            write(*,*)'avant residuals'

            call residualsc(b,npm,ppi,resid_m,pred_m_g,resid_ss &
          ,pred_ss_g,pred_RE)
            ig=0
            ij=0
            do i=1,ns
               do g=1,ng0
                  ig=ig+1
                  ppi0(ig)=PPI(i,g)
               end do
            end do


!            write(*,*)'avant transfo'

            call transfo_estimee(b,npm,nsim0,marker,transfY)

!         else
!            ig=0
!            ij=0
!            do i=1,ns
!               do g=1,ng0
!                  ig=ig+1
!                  ppi0(ig)=0.d0
!               end do
!            end do

         end if

      end if


!      write(*,*)'avant deallocate'

 1589 continue

      deallocate(Y,X,idprob,idea,idg,nmes,prior)

      
      deallocate(zitr,mm,mm1,mm2,im,im1,im2)

      return
      end subroutine hetmixCont



!-----------------------------------------------------------
!                        FUNCPA
!------------------------------------------------------------


      double precision function funcpac(b,npm,id,thi,jd,thj)

      use communc
      use optimc
      IMPLICIT NONE
      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,it,nmestot,ll,ii
      integer ::ier,nmoins,kk
      double precision,dimension(maxmes,nv) ::Z,P,X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: vrais,eps,det,som,thi,thj,temp
      double precision ::Y4,expo,jacobien,beta_densite,ytemp
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3
      double precision,dimension(ng) :: pi
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1

      b1=0.d0
      eps=1.D-20
      do k=1,npm
         b1(k)=b(k)
      end do

      if (id.ne.0) b1(id)=b1(id)+thi
      if (jd.ne.0) b1(jd)=b1(jd)+thj


!----------- rappel des parametres utilises ---------


      Ut=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=b1(nef+j)
               else
                  Ut(j,k)=0.d0
               end if
            end do
         end do
      end if

      If (idiag.eq.0) then
         do j=1,nea
            do k=1,j
               Ut(j,k)=b1(nef+k+j*(+j-1)/2)
            end do
         end do
      end if

! ----------- boucle sur les individus -------------
      nmestot=0
      vrais=0.d0
      jacobien=0.d0


      it=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(it+j,k))
               end do
            end if
         end do
! creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0

         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Se(j,j)=1 ! erreur de mesure = parm de transfo
               Y1(j)=(dble(Y(nmestot))-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))

               jacobien = jacobien - log(b1(nef+nvc+nwg+2))
            end do

         else if (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1

            cc1=abs(b1(nef+nvc+nwg+3))

            dd1=abs(b1(npm))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                  funcpac=-1.d9
                  goto 654
               end if

               jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
               jacobien=jacobien-log(abs(maxY-minY+2*epsY))
            end do

         else if (idlink.eq.2) then ! Splines link



            bb=b1(nef+nvc+nwg+1)

            do kk=2,ntrtot
            splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and. &
                      (Y(nmestot).lt.zitr(kk)))then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                funcpac=-1.d9
                goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+ som +splaa(ll-2)*im2(nmestot)  &
                +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)

               jacobien = jacobien + log(splaa(ll-2)*mm2(nmestot) &
                 +splaa(ll-1)*mm1(nmestot)+splaa(ll)*mm(nmestot))
            end do
         end if


!         if (i.lt.3)then
!            write(*,*)'nmes',nmes(i),b1((nef+nvc+nwg+1):npm),nef
!            write(*,*)'Y1',Y1
!            write(*,*)'Y',Y(nmestot-nmes(i)+1:nmestot)
!         end if


! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Se

! Vi en vecteur

            jj=0
            Vi=0.d0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL dsinvc(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               funcpac=-1.d9
               goto 654
            end if

!     retransformation du vecteur Vi en matrice :

            VC=0.d0
            do j=1,nmes(i)
               do k=1,nmes(i)
                  if (k.ge.j) then
                     VC(j,k)=Vi(j+k*(k-1)/2)
                  else
                     VC(j,k)=Vi(k+j*(j-1)/2)
                  end if
               end do
            end do
         end if
!     debut du calcul de la vraisemblance
       vrais=vrais-nmes(i)*dlog(dble(2*3.14159265))
! contribution individuelle a la vraisemblance

! cas 1 : ng=1

       if (ng.eq.1) then

            vrais=vrais-det
            b0=0.d0
            l=0
            X00=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
		     X00(j,l)=dble(X(it+j,k))
                  end do
! idg ne 0 pour l'iintercept forcement donc on met le parm a 0
                  if (k.eq.1) then
                     b0(l)=0.d0
                  else
                     b0(l)=b1(nprob+l-1)
                  end if
                end if
            end do


            mu=0.d0
            y2=0.d0
            y3=0.d0
            y4=0.d0
            mu=matmul(X00,b0)
            Y2=Y1-mu
            Y3=matmul(VC,Y2)
            Y4=DOT_PRODUCT(Y2,Y3)

            vrais=vrais-Y4

! cas 2 :  ng>1  composantes
         else

            if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
            else

! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
                Xprob=0.d0
                Xprob(1)=1
                l=0
                do k=1,nv
                    if (idprob(k).eq.1) then
                    l=l+1
                    Xprob(1+l)=X(it+1,k)
                    end if
                end do
                pi=0.d0
                temp=0.d0
                Do g=1,ng-1
                    bprob=0.d0
                    do k=1,nvarprob
                        bprob(k)=b1((k-1)*(ng-1)+g)
                    end do

                temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

            end do

            pi(ng)=1/(1+temp)

            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do

          end if

! creation des vecteurs de variables explicatives
            l=0
            m=0
            X00=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
			X2(j,l)=dble(X(it+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X00(j,m)=dble(X(it+j,k))
                  end do
               end if
            end do

            b2=0.d0
            b0=0.d0
            expo=0.d0
            do g=1,ng
               nmoins=0
               l2=0
               m2=0
               do k=1,nv
                  if (idg(k).eq.1) then
! parametre a 0 pour l'intercept
                     if (k.eq.1) then
                        m2=m2+1
                        b0(m2)=0.d0
                     else
                        m2=m2+1
                        b0(m2)=b1(nprob+nmoins+1)
                        nmoins=nmoins+1
                     end if
                  else if (idg(k).eq.2) then
! parametre a 0 pour l'intercept de la premiere classe
                     if (k.eq.1) then
                        if (g.eq.1) then
                            l2=l2+1
                            b2(l2)=0.d0
                            nmoins=nmoins+ng-1
                        else
                            l2=l2+1
                            b2(l2)=b1(nprob+nmoins+g-1)
                            nmoins=nmoins+ng-1
                        end if
                     else
                        l2=l2+1
                        b2(l2)=b1(nprob+nmoins+g)
                        nmoins=nmoins+ng
                     end if
                  end if
               end do


!               if (thi.eq.0.and.thj.eq.0) then
!                if (i.eq.1) then
!                    write(*,*)'g',g,b2
!                    write(*,*)'g',g,b0
!                    stop
!                    end if
!                    end if


! variance covariance si spec aux classes :

               Ut1=Ut
               if (nwg.ne.0) then
                  Ut1=0.d0
                  if (g.eq.ng) then
                     Ut1=Ut
                  else
                     Ut1=Ut*abs(b1(nef+nvc+g))
                  end if


                  P=0.d0
                  P=MATMUL(Z,Ut1)
                  VC=0.d0
                  VC=MATMUL(P,transpose(P))+Se

! Vi en vecteur
                  Vi=0.d0
                  jj=0
                  do j=1,nmes(i)
                     do k=j,nmes(i)
                        jj=j+k*(k-1)/2
                        Vi(jj)=VC(j,k)
                     end do
                  end do

                  CALL dsinvc(Vi,nmes(i),eps,ier,det)
                  if (ier.eq.-1) then
                     funcpac=-1.d9
                     goto 654
                  end if

!     retransformation du vecteur Vi en matrice :

                  VC=0.d0
                  do j=1,nmes(i)
                     do k=1,nmes(i)
                        if (k.ge.j) then
                           VC(j,k)=Vi(j+k*(k-1)/2)
                        else
                           VC(j,k)=Vi(k+j*(j-1)/2)
                        end if
                     end do
                  end do

               end if

               mu=0.d0
               y2=0.d0
               y3=0.d0
               y4=0.d0
               mu=matmul(X00,b0)+matmul(X2,b2)
               Y2=Y1-mu
               Y3=Matmul(VC,Y2)
               Y4=0.d0
               Y4=DOT_PRODUCT(Y2,Y3)

               expo = expo+pi(g)*exp((-det-Y4)/2.d0)

            end do
            vrais=vrais+2*log(expo)

         end if

         it=it+nmes(i)
       end do

! FIN BOUCLE SUJET

      funcpac=vrais/2.D0+jacobien

 654  continue

      return

      end function funcpac


!------------------------------------------------------------
!                      POSTPROB
!------------------------------------------------------------

!-------------------------------------------------------------
!
!          Subroutine pour calculer les
!      probabilites a posteriori de suivre chacune
!      des composantes g pour chacun des sujets i
!
!-------------------------------------------------------------

      subroutine postprobc(b,npm,PPI)
      use communc
      use optimc
      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk,nmestot,ii,ll
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det,som,temp,Y4,f,ytemp
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1


      eps=1.D-20

      PPI=0.D0

      do k=1,npm
         b1(k)=b(k)
      end do


!----------- rappel des parametres utilises ---------


      Ut=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=b1(nef+j)
               else
                  Ut(j,k)=0.d0
               end if
            end do
         end do
      end if

      If (idiag.eq.0) then
         do j=1,nea
            do k=1,j
               Ut(j,k)=b1(nef+k+j*(j-1)/2)
            end do
         end do
      end if


! ----------- boucle sur les individus -------------
      nmestot=0
      it=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+se*seIni ----------

! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(it+j,k))
               end do
            end if
         end do



! creation de s2*I et Y1
         Se=0.d0
         Y1=0.d0
         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Se(j,j)=1.d0
               Y1(j)=dble(Y(nmestot)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
            end do

         else if (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
              (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
              (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                PPI=-1.d0
                go to 147
               end if

            end do

         else if (idlink.eq.2) then ! Splines link

            bb=b1(nef+nvc+nwg+1)


            do kk=2,ntrtot
               splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and.(Y(nmestot).lt.zitr(kk)))then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                PPI=-1.d0
                go to 147
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+som +splaa(ll-2)*im2(nmestot) &
               +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)
            end do
        end if

! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

        if (nwg.eq.0) then
          P=0.d0
          P=MATMUL(Z,Ut)
          VC=0.d0
          VC=MATMUL(P,transpose(P))+Se

! Vi en vecteur
         jj=0
         do j=1,nmes(i)
            do k=j,nmes(i)
               jj=j+k*(k-1)/2
               Vi(jj)=VC(j,k)
            end do
         end do

         CALL dsinvc(Vi,nmes(i),eps,ier,det)
         if (ier.eq.-1) then
            PPI=-1.d0
            go to 147
         end if

! retransformation du vecteur Vi en matrice :

         do j=1,nmes(i)
            do k=1,nmes(i)
               if (k.ge.j) then
                  VC(j,k)=Vi(j+k*(k-1)/2)
               else
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end if
            end do
         end do
       end if

       if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
       else

! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       Xprob=0.d0
       Xprob(1)=1
       l=0
       do k=1,nv
          if (idprob(k).eq.1) then
             l=l+1
             Xprob(1+l)=X(it+1,k)
          end if
       end do

       pi=0.d0
       temp=0.d0
       Do g=1,ng-1
          bprob=0.d0
          do k=1,nvarprob
             bprob(k)=b1((k-1)*(ng-1)+g)
          end do

          temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

          pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

       end do

       pi(ng)=1/(1+temp)

       do g=1,ng-1
          pi(g)=pi(g)*pi(ng)
       end do

       end if
!     creation des vecteurs de variables explicatives
       l=0
       m=0
       X0=0.d0
       X2=0.d0
       do k=1,nv
          if (idg(k).eq.2) then
             l=l+1
             do j=1,nmes(i)
		X2(j,l)=dble(X(it+j,k))
             end do
          else if (idg(k).eq.1) then
             m=m+1
             do j=1,nmes(i)
                X0(j,m)=dble(X(it+j,k))
             end do
          end if
       end do
!     calcul de la vraisemblance par composante
       f=0.d0
       fi=0.d0
       do g=1,ng
          nmoins=0
          l2=0
          m2=0
          do k=1,nv
             if (idg(k).eq.1) then
                ! parametre a 0 pour l'intercept
                if (k.eq.1) then
                   m2=m2+1
                   b0(m2)=0.d0
                else
                   m2=m2+1
                   b0(m2)=b1(nprob+nmoins+1)
                   nmoins=nmoins+1
                end if
             else if (idg(k).eq.2) then
                ! parametre a 0 pour l'intercept de la premiere classe
                 if (k.eq.1) then
                        if (g.eq.1) then
                            l2=l2+1
                            b2(l2)=0.d0
                            nmoins=nmoins+ng-1
                        else
                            l2=l2+1
                            b2(l2)=b1(nprob+nmoins+g-1)
                            nmoins=nmoins+ng-1
                        end if
                else
                   l2=l2+1
                   b2(l2)=b1(nprob+nmoins+g)
                   nmoins=nmoins+ng
                end if
             end if
          end do
!     variance covariance si spec aux classes :
          if (nwg.ne.0) then
             Ut1=0.d0
             if (g.eq.ng) then
                Ut1=Ut
             else
                Ut1=Ut*abs(b1(nef+nvc+g))
             end if

             P=0.d0
             P=MATMUL(Z,Ut1)
             VC=0.d0
             VC=MATMUL(P,transpose(P))+Se
!     Vi en vecteur
             jj=0
             do j=1,nmes(i)
                do k=j,nmes(i)
                   jj=j+k*(k-1)/2
                   Vi(jj)=VC(j,k)
                end do
             end do

             CALL dsinvc(Vi,nmes(i),eps,ier,det)
             if (ier.eq.-1) then
                PPI=-1.d0
                goto 147
             end if

!     retransformation du vecteur Vi en matrice :

             do j=1,nmes(i)
                do k=1,nmes(i)
                   if (k.ge.j) then
                      VC(j,k)=Vi(j+k*(k-1)/2)
                   else
                      VC(j,k)=Vi(k+j*(j-1)/2)
                   end if
                end do
             end do

          end if

          mu=0.d0
          mu=matmul(X0,b0)+matmul(X2,b2)

          Y2=Y1-mu
          Y3=MATMUL(VC,Y2)
          Y4=DOT_PRODUCT(Y2,Y3)
          fi(g)=fi(g)- nmes(i)*log(dble(2*3.14159265))
          fi(g)=fi(g) -det
          fi(g)=fi(g) - Y4
          fi(g)=fi(g)/(2.d0)
          fi(g)=exp(fi(g))

       end do
       f=DOT_PRODUCT(pi,fi)
       do g=1,ng
          PPI(i,g)=pi(g)*fi(g)/f
       end do

       it=it+nmes(i)

      end do

 147  continue
      return

      end subroutine postprobc

!------------------------------------------------------------
!                      RESIDUALS
!------------------------------------------------------------


      subroutine residualsc(b1,npm,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE)

      use communc

      use optimc

      implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,npm,nmestot,ll,ii
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea) ::err2
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se,VC1
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob,b3
      double precision :: eps,det,temp
      double precision,dimension(nea,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m,resid_ss
      double precision,dimension(nobs*ng)::pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(ns*nea)::pred_RE
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,ytemp,som,cc1

      eps=1.D-20

!----------- rappel des parametres utilises ---------

! creation de Ut, decomposition de cholesky pour G
      Ut=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=b1(nef+j)
               else
                  Ut(j,k)=0.d0
               end if
            end do
         end do
      end if

      If (idiag.eq.0) then
         do j=1,nea
            do k=1,j
               Ut(j,k)=b1(nef+k+j*(j-1)/2)
            end do
         end do
      end if

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m_g=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0

      nmes_cur=0
      kk=0
      nmestot=0
      pred_RE=0.d0
      do i=1,ns

!     -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(nmes_cur+j,k))
               end do
            end if

         end do
! creation de s2*I et Y1
         Se=0.d0
         Y1=0.d0
         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Se(j,j)=1.d0
               Y1(j)=dble(Y(nmestot)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
            end do

         elseif (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1


            do j=1,nmes(i)

               nmestot=nmestot+1

               Se(j,j)=1.d0
               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1

               if (Y1(j).eq.999.d0) then
                  do k=1,nmes(i)
                     resid_m(nmes_cur+k)=9999.d0
                     pred_m_g(nmes_cur+k)=9999.d0

                     resid_ss(nmes_cur+k)=9999.d0
                     pred_ss_g(nmes_cur+k)=9999.d0
                  end do
                  do k=1,nea
                     pred_RE((i-1)*nea+k)=9999.d0
                  end do
                  goto 654
               end if

            end do

         elseif (idlink.eq.2) then ! Splines link


            bb=b1(nef+nvc+nwg+1)

            do kk=2,ntrtot
               splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and.(Y(nmestot).lt.zitr(kk)))then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
               do k=1,nmes(i)
                  resid_m(nmes_cur+k)=9999.d0
                  pred_m_g(nmes_cur+k)=9999.d0

                  resid_ss(nmes_cur+k)=9999.d0
                  pred_ss_g(nmes_cur+k)=9999.d0
               end do
               do k=1,nea
                    pred_RE((i-1)*nea+k)=9999.d0
                end do
               goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+som +splaa(ll-2)*im2(nmestot) &
                 +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)

            end do

         end if
!     debut du calcul de la vraisemblance

!     cas 1 : ng=1

         if (ng.eq.1) then


            Valea=0.d0
            VC=0.d0
            P=0.d0


            P=MATMUL(Z,Ut)
            Valea=MATMUL(Ut,transpose(P))
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Se

!     Vi en vecteur

            jj=0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL dsinvc(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               do j=1,nmes(i)
                  resid_m(nmes_cur+j)=9999.d0
                  pred_m_g(nmes_cur+j)=9999.d0

                  resid_ss(nmes_cur+j)=9999.d0
                  pred_ss_g(nmes_cur+j)=9999.d0
               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
               end do
               goto 654
            end if

!     retransformation du vecteur Vi en matrice :
            VC1=0.d0
            do j=1,nmes(i)
               do k=1,nmes(i)
                  if (k.ge.j) then
                     VC1(j,k)=Vi(j+k*(k-1)/2)
                  else
                     VC1(j,k)=Vi(k+j*(j-1)/2)
                  end if
               end do
            end do


            b0=0.d0
            l=0
            X0=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
                     X0(j,l)=dble(X(nmes_cur+j,k))
                  end do
                  if (k.gt.1) then
                     b0(l)=b1(nprob+l-1)
                  end if
               end if
            end do

            mu=matmul(X0,b0)
            Y2=Y1-mu

            err1=0.d0
            err1=MATMUL(VC1,Y2)
            err2=0.d0
            err2=MATMUL(Valea,err1)
            pred1=0.d0
            pred1=mu+MATMUL(Z,err2)
            err1=0.d0
            err1=MATMUL(Z,err2)
            do j=1,nmes(i)
               resid_m(nmes_cur+j)=Y2(j)
               pred_m_g(nmes_cur+j)=mu(j)

               resid_ss(nmes_cur+j)=Y2(j)-err1(j)
               pred_ss_g(nmes_cur+j)=pred1(j)
            end do


            do k=1,nea
               pred_RE((i-1)*nea+k)=err2(k)
            end do


!     cas 2 :  ng>1  composantes
         else

            if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
            else

!     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            Xprob=0.d0
            Xprob(1)=1.d0
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(nmes_cur+1,k)
               end if
            end do
!     write(*,*)'l apres Xprob',l,(Xprob(j),j=1,10)
            pi=0.d0
            temp=0.d0
            Do g=1,ng-1
               bprob=0.d0
               do k=1,nvarprob
                  bprob(k)=b1((k-1)*(ng-1)+g)
               end do
               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

            end do
            pi(ng)=1.d0/(1.d0+temp)
            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do


           end if
!     creation des vecteurs de variables explicatives
            l=0
            m=0
            X0=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
                     X2(j,l)=dble(X(nmes_cur+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X0(j,m)=dble(X(nmes_cur+j,k))
                  end do
               end if
            end do


            do g=1,ng
               nmoins=0
               l2=0
               m2=0
               n2=0
               b0=0.d0
               b2=0.d0
               b3=0.d0
               nmoins2=0

               do k=1,nv
                  if (idg(k).eq.1) then
! parametre a 0 pour l'intercept
                     if (k.eq.1) then
                        m2=m2+1
                        b0(m2)=0.d0
                     else
                        m2=m2+1
                        b0(m2)=b1(nprob+nmoins+1)
                        nmoins=nmoins+1
                     end if
                  else if (idg(k).eq.2) then
                     if (k.eq.1) then
                        if (g.eq.1) then
                            l2=l2+1
                            b2(l2)=0.d0
                            nmoins=nmoins+ng-1
                        else
                            l2=l2+1
                            b2(l2)=b1(nprob+nmoins+g-1)
                            nmoins=nmoins+ng-1
                        end if
                     else
                        l2=l2+1
                        b2(l2)=b1(nprob+nmoins+g)
                        nmoins=nmoins+ng
                     end if
                  end if
                  IF (IDEA(k).EQ.1.and.idg(k).eq.1) THEN
                  ! parametre a 0 pour l'intercept
                     if (k.eq.1) then
                        n2=n2+1
                        b3(n2)=0.d0
                     else
                        n2=n2+1
                        b3(n2)=b1(nprob+nmoins2+1)
                        nmoins2=nmoins2+1
                     end if
                  else if(IDEA(k).EQ.1.and.idg(k).eq.2) THEN
                    if (k.eq.1) then
                        if (g.eq.1) then
                            n2=n2+1
                            b3(n2)=0.d0
                            nmoins2=nmoins2+ng-1
                        else
                            n2=n2+1
                            b3(n2)=b1(nprob+nmoins2+g-1)
                            nmoins2=nmoins2+ng-1
                        end if
                     else
                        n2=n2+1
                        b3(n2)=b1(nprob+nmoins2+g)
                        nmoins2=nmoins2+ng
                     end if
                  end if
               end do

               mu=0.d0
               mu=matmul(X0,b0)+matmul(X2,b2)

               VC=0.d0
               P=0.d0
               Ut1=Ut
               if (nwg.ne.0) then
                  if (g.eq.ng) then
                     Ut1=Ut
                  else
                     Ut1=Ut*b1(nef+nvc+g)
                  end if
               end if
               P=0.d0
               Valea=0.d0
               VC=0.d0
               P=MATMUL(Z,Ut1)
               Valea=MATMUL(Ut1,transpose(P))
               VC=MATMUL(P,transpose(P))+Se
!     Vi en vecteur
               jj=0
               do j=1,nmes(i)
                  do k=j,nmes(i)
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do
               CALL dsinvc(Vi,nmes(i),eps,ier,det)
               if (ier.eq.-1) then
                  do j=1,nmes(i)
                     resid_m(nmes_cur+j)=9999.d0
                     resid_ss(nmes_cur+j)=9999.d0
                     do l=1,ng
                        pred_ss_g((l-1)*nobs+nmes_cur+j)=9999.d0
                        pred_m_g((l-1)*nobs+nmes_cur+j)=9999.d0
                     end do
                  end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
               end do
                  goto 654
               end if

!     retransformation du vecteur Vi en matrice :
               VC1=0.d0
               do j=1,nmes(i)
                  do k=1,nmes(i)
                     if (k.ge.j) then
                        VC1(j,k)=Vi(j+k*(k-1)/2)
                     else
                        VC1(j,k)=Vi(k+j*(j-1)/2)
                     end if
                  end do
               end do


               Y2=Y1-mu

               err1=0.d0
               err1=MATMUL(VC1,Y2)
               err2=0.d0
               err2=MATMUL(Valea,err1)
               pred1=0.d0
               pred1=mu+MATMUL(Z,err2)

               do j=1,nmes(i)
                  pred_m_g((g-1)*nobs+nmes_cur+j)=mu(j)
                  pred_ss_g((g-1)*nobs+nmes_cur+j)=pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j) &
                      +ppi(i,g)*(Y1(j)-pred1(j))
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))

               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*err2(k)
               end do

            end do

         end if

 654  continue

         nmes_cur=nmes_cur+nmes(i)
      end do

! FIN BOUCLE SUJET


      end subroutine residualsc





! =============================================
! subroutine de creation de design matrix
! =============================================




      subroutine design_splines (ier)

      use communc

      implicit none
      
      integer ::j,i,jj,l,k,ier
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

      ier=0
      jj=0
      l=0
      do i=1,ns
         do j=1,nmes(i)
            jj=jj+1
!     ou se trouve la valeur de zi

            do k = 2,ntrtot-2
               if ((Y(jj).ge.zitr(k-1)).and.(Y(jj).lt.zitr(k)))then
                  l=k-1
               endif

            end do


            if (Y(jj).eq.zitr(ntrtot-2)) then
               l=ntrtot-3
            end if

            ht2 = zitr(l+1)-Y(jj)
            htm= Y(jj)-zitr(l-1)
            ht = Y(jj)-zitr(l)
            ht3 = zitr(l+2)-Y(jj)
            hht = Y(jj)-zitr(l-2)
            h = zitr(l+1)-zitr(l)
            hh= zitr(l+1)-zitr(l-1)
            hn= zitr(l+1)-zitr(l-2)
            h2n=zitr(l+2)-zitr(l-1)
            h2= zitr(l+2)-zitr(l)
            h3= zitr(l+3)-zitr(l)

            if (Y(jj).ne.zitr(ntrtot-2)) then
               mm2(jj) = (3.d0*ht2*ht2)/(hh*h*hn)
               mm1(jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
               mm(jj)  = (3.d0*ht*ht)/(h3*h2*h)

            end if
            if (Y(jj).eq.zitr(ntrtot-2)) then
               mm2(jj) = 0.d0
               mm1(jj) = 0.d0
               mm(jj)  = 3.d0/h
            end if

            if (mm2(jj).lt.0.or.mm1(jj).lt.0.or.mm(jj).lt.0)then
                ier=-1
                goto 765
            end if

            im2(jj)=hht*mm2(jj)/(3.d0)+ h2n*mm1(jj)/(3.d0) &
             +h3*mm(jj)/(3.d0)
            im1(jj)=htm*mm1(jj)/(3.d0)+h3*mm(jj)/(3.d0)
            im(jj)=ht*mm(jj)/(3.d0)

         end do
      end do


765     continue

      end subroutine design_splines






! =================================================================
!              Densite d'une beta
!=================================================================

      double precision Function beta_densite(X,a,b)

      implicit none
      
      double precision :: beta,a,b,gammln,X
 
      beta=exp(gammln(a+b)-gammln(a)-gammln(b))
      
      beta_densite=((X)**(a-1))*((1-X)**(b-1))*beta

      return
      
      end Function beta_densite


! Calcul de gamma (gammln(a))

      double precision Function gammln(xx)

! retourne la valeur ln(gamma(xx)) pour xx>0
      implicit none 
      
      integer::j
      double precision:: ser,stp,tmp,x,y,xx
      double precision,dimension(6)::cof
      
      save cof,stp

      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
      24.01409824083091d0, -1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      end do
      gammln=tmp+log(stp*ser/x)

      return
      
      end Function gammln






! =================================================================
!              CDF incomplete d'une beta
!=================================================================


      double precision Function betai(a,b,x)

      implicit none
      
      double precision :: a,b,x,bt,betaCF,gammln,temp
      
      if(x.lt.0.d0.or.x.gt.1.d0) then
            betai=999.d0
            return
      end if
      if (x.eq.0..or.X.eq.1.) then
         bt=0.
      else
         bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      end if

      if (x.lt.(a+1.)/(a+B+2.)) then
         temp=betaCF(a,b,x)
         if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=bt*temp/A
         return
      else
        temp=betaCF(b,a,1.-x)
        if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=1.-bt*temp/b
         return
      end if
      
      end Function betai

! betaCF est utilis� par betai

      double precision Function betaCF(a,b,x)

      implicit none
      
      integer ::m,m2
      integer,parameter ::maxit=100
      double precision ::a,b,x,aa,c,d,del,h,qab,qam,qap
      double precision,parameter::eps=3.e-7,fpmin=1.e-30

      
      qab=a+b
      qap=a+1
      qam=a-1
      c=1.       ! first step of Lentz's method
      d=1.-qab*x/qap
      if (abs(d).lt.fpmin) d=fpmin
      d=1./d
      h=d
      do m=1,maxit
         m2=2*m
         aa=m*(b-m)*x/((qam+m2)*(a+m2))
         d=1.+aa*d        ! one step (the even one) of the recurrence
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         h=h*d*c
         aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
         d=1.+aa*d        ! next step of the recurrence (the odd one)
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h=h*del
         if (abs(del-1.).lt.eps) goto 1
      end do

!  pause 'a or b too big, or maxit too small in betaCF'
            betaCF=999.d0
!            write(*,*)'problem'
            return

 1    betaCF=h
      return
      end Function betaCF




!---------------------------------------------------------------

      !transfo estimee

!---------------------------------------------------------------

      subroutine transfo_estimee(b,npm,nsim,marker,transfY)

      use communc
      
      implicit none
      
      double precision,dimension(nsim)::marker,transfY
      double precision,dimension(ntrtot)::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,eps,pas,ytemp,cc1
      double precision, dimension(npm)::b,b1
      integer::kk,nsim,npm,j,k



      b1=0.d0
      eps=1.D-20
      do k=1,npm
         b1(k)=b(k)
      end do


!       write(*,*)'infos',minY,maxY,nsim,npm
!       write(*,*)'b',(b1(j),j=1,npm)



       marker=0.d0
       transfY=0.d0

       pas=(maxY-minY)/dble(nsim-1)
       j=1
       marker(1)=minY
       do while(j.lt.nsim)
           j=j+1
           marker(j)=marker(j-1)+pas
       end do
       marker(nsim)=maxY

!       write(*,*)(marker(j),j=1,nsim)

       if (idlink.eq.2) then
        splaa=0.d0


        splaa(1)=b1(nef+nvc+nwg+1)
        do kk=2,ntrtot
           splaa(kk)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
             end do

!           write(*,*)'ntr',(b1(nprob+nrisq
!     &       +nvarxevt+ncssg+ncg*ng+kk),kk=1,ntrtot)


            call   estim_splines_ssstd(nsim,splaa,marker,transfY)

        else if (idlink.eq.1) then

            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nsim
                  ytemp=(marker(j)-minY+epsY)/(maxY-minY+2*epsY)
                  transfY(j)=(betai(aa,bb,ytemp)-cc1)/dd1
                  if (transfY(j).eq.999.d0) then
!                    write(*,*)'problem'
                  end if

               end do


        else if (idlink.eq.0) then

                 do j=1,nsim
                    transfY(j)=(marker(j)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
                 end do
        end if
        end subroutine transfo_estimee





!---------------------------------------------------------------

!      TRANSFOS SPLINES SANS STD ERR

!---------------------------------------------------------------





        subroutine estim_splines_ssstd(nsim,aa,test,transf)

       use communc
       
       implicit none
       
       integer::nsim,j,k,l
       double precision,dimension(nsim)::mmm,mmm1,mmm2 &
      ,iim,iim1,iim2,transf,test
       double precision,dimension(ntrtot)::aa,Xspl
       double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn


! matrice de transition pour delta-metho (carre des parms 2,..,ntr)

       do j=1,nsim
! ou se trouve la valeur
         l=0

         do k = 2,ntrtot-2
               if ((test(j).ge.zitr(k-1)).and.(test(j).lt.zitr(k)))then
                  l=k-1
               endif
            end do

           if (test(j).eq.zitr(ntrtot-2)) then
               l=ntrtot-3
            end if

!         if (l.lt.1.or.l.gt.ntrtot-1) then
!            write(*,*)'probleme estim splines',l
!            write(*,*)'j=',j,'test(j)',test(j)
!            stop
!         end if


               ht2 = zitr(l+1)-test(j)
               htm= test(j)-zitr(l-1)
               ht = test(j)-zitr(l)
               ht3 = zitr(l+2)-test(j)
               hht = test(j)-zitr(l-2)
               h = zitr(l+1)-zitr(l)
               hh= zitr(l+1)-zitr(l-1)
               hn= zitr(l+1)-zitr(l-2)
               h2n=zitr(l+2)-zitr(l-1)
               h2= zitr(l+2)-zitr(l)
               h3= zitr(l+3)-zitr(l)

               if (test(j).ne.zitr(ntrtot-2)) then
                  mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
                  mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
               end if
               if (test(j).eq.zitr(ntrtot-2)) then
                  mmm2(j) = 0.d0
                  mmm1(j) = 0.d0
                  mmm(j)  = 3.d0/h
               end if

               iim2(j)=hht*mmm2(j)/(3.d0)+ h2n*mmm1(j)/(3.d0) &
      +h3*mmm(j)/(3.d0)
               iim1(j)=htm*mmm1(j)/(3.d0)+h3*mmm(j)/(3.d0)

               iim(j)=ht*mmm(j)/(3.d0)

!-------- transformation et IC de la transformation :

		    Xspl=0.d0
            Xspl(1)=1
            do k=2,l
               Xspl(k)=1
            end do
            Xspl(l+1)=iim2(j)
            Xspl(l+2)=iim1(j)
            Xspl(l+3)=iim(j)
            transf(j)= dot_product(Xspl,aa)
      end do



      end subroutine estim_splines_ssstd

