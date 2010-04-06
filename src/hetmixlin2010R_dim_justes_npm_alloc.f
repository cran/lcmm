c==========================================================
c
c      Program for estimating linear mixed models
c         with a mixture of distribution for
c                  the random effects
c
c        Cecile Proust, Helene Jacqmin-Gadda
c
c
c       Corresponding author :
c       Cécile Proust, INSERM U897, ISPED,
c       146  rue L\'eo Saignat,
c       33076 Bordeaux cedex, France.
c       Tel: (33) 5 57 57 45 79; Fax: (33) 5 56 24 00 81;
c       e-mail : cecile.proust@isped.u-bordeaux2.fr
c
c                                        10/12/08
c===========================================================
c
c
c--------- MAIN VARIABLES  ---------------------------------
c
c nv : number of covariates in the model
c ng : number of components in the model
c ns : number of units
c Numero(i) : identification number of unit i
c Y(i,j) : dependent variable for unit i at time j
c X(i,j,k) : covariate k for unit i at time j
c idea(k) : indicator of random effect on covariate k :
c                =0 if no random effect
c                =1 if random effect
c idg(k) : indicator of mixture on covariate k :
c                =0 if no effect
C                =1 if overall effect
c                =2 if class-specific effect
c idprob(k) : indicator that covariate k is in the class-membership model :
c                =0 if no effect
C                =1 if overall effect
c                =2 if class-specific effect
c idiag : indicator of the random-effect covariance matrix
c         structure
c                =0 if unstructured matrix
c                =1 if diagonal matrix
c
c
c --------- REMARKS ---------------------------------------
c
c Transformation for the covariance matrix V : Cholesky matrix U
c     => V=U'U
c
c multinomial logistic model for class probabilities
c
c Algorithm of optimisation : modified Marquardt
c (if the loglikelihood at iteration k is not improved
c then the step for the update is modified to ensure
c that the loglikelihood is improved)
c
c
c----------- PROGRAM COMPOSITION --------------------------
c
c     - SUBOURTINE HETMIXLIN
c     - function FUNCPA
c     - subroutine MARQ9
c     - subroutine DERIVA
c     - subroutine SEARPAS
c     - subroutine VALFPA
c     - function MAXT
c     - subroutine DCHOLE
c     - subroutine DMFSD
c     - subroutine DSINV
c     - subroutine POSTPROB
c
c----------------------------------------------------------


      module parameters

      double precision,parameter::epsa=1.d-3,epsb=1.d-3,epsd=1.d-3
      end module parameters

      module commun

      use parameters

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg
     & ,nprob,nvarprob,maxmes,nobs
      double precision,dimension(:,:),allocatable,save::Y
      double precision,dimension(:,:,:),allocatable,save ::X
      integer,dimension(:),allocatable,save ::idea,idg,idprob
      integer,dimension(:),allocatable,save :: nmes

      end module commun



C================ SUBROUTINES ================================

      subroutine hetmixlin(
     & Y0,X0,
     &     idprob0,idea0,idg0
     &     ,ns0,ng0,nv0,nobs0
     &     ,nmes0,idiag0,nwg0
     &     ,npm0,b
     & ,se
     &     ,vrais
     &     ,ni,istop,gconv
     &,ppi0
     &     ,resid_m0,resid_ss0
     &      ,pred_m_g0
     &     ,pred_ss_g0
     & )


      use commun
      implicit none


      integer::ns0,ng0,nv0,nobs0,idiag0,nwg0,jtemp,npm0
      double precision,dimension(npm0,npm0) ::VC,Ut,U
      double precision,dimension(ns0*ng0) ::ppi0
      double precision,dimension(ns0,ng0) ::PPI
      integer :: nef,i,g,j,ij,npm,ier,istop,id,ni,k,ktemp,ig
     &     ,nmestot
      double precision,dimension(npm0)::b,se,mvc
      double precision ::vrais,eps
      double precision,dimension(npm0*(npm0+3)/2) ::V
      double precision,dimension(3)::gconv

      double precision,dimension(nobs0)::Y0
      double precision,dimension(nobs0*nv0)::X0
      integer,dimension(nv0)::idea0,idg0,idprob0
      integer,dimension(ns0):: nmes0
      double precision :: ca,cb,dd
      double precision,dimension(nobs0)::resid_m0
     &     ,resid_ss0
      double precision,dimension(nobs0)::resid_m
     &     ,resid_ss
      double precision,dimension(nobs0*ng0)::
     &     pred_m_g0,pred_ss_g0
      double precision,dimension(nobs0,ng0)::
     &     pred_m_g,pred_ss_g


      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do


      allocate(Y(ns0,maxmes),X(ns0,maxmes,nv0),idprob(nv0)
     & ,idea(nv0),idg(nv0),nmes(ns0))


      ppi0=0.d0
      resid_m0=0.d0
      resid_ss0=0.d0
      pred_m_g0=0.d0
      pred_ss_g0=0.d0
      se=0.d0

      eps=1.d-20
C enrigstrement pour les modules
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
         DO i=1,ns
            if (k.eq.1) then
               nmes(i)=nmes0(i)
               do j=1,nmes(i)
                  nmestot=nmestot+1
                  jtemp=jtemp+1
                  Y(i,j)=Y0(jtemp)
               end do
            end if
            do j=1,nmes(i)
               ktemp=ktemp+1
               X(i,j,k)=X0(ktemp)
            end do
         end do
      end do

C creation des parametres

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


c  nb effets fixes = nb effets fixes sans melange
c                  + ng fois le nb de var dans melange


      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

      nef=nprob+ncssg+ncg*ng
      npm=nef+nvc+nwg+1


      if (idiag.eq.1) then
         DO j=1,nvc
            B(nef+j)=dsqrt(abs(B(nef+j)))
         END DO
      end if

c si idiag=0, on met dans le vecteur des parms, les parms
c de la transformee de Cholesky

      if (idiag.eq.0) then


         DO j=1,nvc
            mvc(j)=B(nef+j)
         END DO

         CALL DMFSD(mvc,nea,EPS,IER)
         DO j=1,nvc
            B(nef+j)=mvc(j)
         END DO
      end if

      B(nef+nvc+nwg+1)=abs(B(nef+nvc+nwg+1))

C lancement de l'optimisation

      IF (npm.eq.1) then
         istop=10
         go to 1589
      else
         ca=0.d0
         cb=0.d0
         dd=0.d0

         call marq98(b,npm,ni,V,vrais,ier,istop,ca,cb,dd)

C         write(*,*)
C         write(*,*)'    FIN OPTIMISATION  ..... '
C         write(*,*)'istop',istop,'vrais',vrais

         gconv=0.d0
         gconv(1)=ca
         gconv(2)=cb
         gconv(3)=dd

         do j=1,npm
            se(j)=0.d0
         end do

         do k=1,nwg+1
            b(nef+nvc+k)=abs(b(nef+nvc+k))
         end do

C resultats : SE parms transformés

         if (istop.eq.1.or.istop.eq.3) then
            id=1
            DO i=1,npm
               IF(istop.eq.1) then
                  se(i)=dsqrt(V(id))
               else
                  se(i)=0.d0
               end if
               id=id+i+1
            end do
            do i=1,nvc
               se(nef+i)=0.d0
            end do

C probas posteriori

            if (ng.gt.1) then
               call postprob(B,npm,PPI)
            end if



            call residuals(b,npm,ppi,resid_m,pred_m_g,resid_ss
     &     ,pred_ss_g)


            ig=0
            ij=0
            do i=1,ns
               do g=1,ng0
                  ig=ig+1
                  ppi0(ig)=PPI(i,g)
               end do
               do j=1,nmes(i)
                  ij=ij+1
                  resid_ss0(ij)=resid_ss(ij)
                  resid_m0(ij)=resid_m(ij)
                  do g=1,ng0
                     pred_ss_g0(ij+nmestot*(g-1))=pred_ss_g(ij,g)
                     pred_m_g0(ij+nmestot*(g-1))=pred_m_g(ij,g)
                  end do
               end do
            end do

         else


            ig=0
            ij=0
            do i=1,ns
               do g=1,ng0
                  ig=ig+1
                  ppi0(ig)=0.d0
               end do
               do j=1,nmes(i)
                  ij=ij+1
                  resid_ss0(ij)=0.d0
                  resid_m0(ij)=0.d0
                  do g=1,ng0
                     pred_ss_g0(ij+nmestot*(g-1))=0.d0
                     pred_m_g0(ij+nmestot*(g-1))=0.d0
                  end do
               end do
            end do

         end if

         U=0.d0
         Ut=0.d0
         If (idiag.eq.1) then
            do j=1,nea
               U(j,j)=b(nef+j)
               Ut(j,j)=b(nef+j)
            end do
         end if
         If (idiag.eq.0) then
            do j=1,nea
               do k=j,nea
                  U(j,k)=b(nef+j+k*(k-1)/2)
                  Ut(k,j)=b(nef+j+k*(k-1)/2)
               end do
            end do
         end if


         VC=Matmul(Ut,U)


         do j=1,nea
            If (idiag.eq.1) then
               b(nef+j)=VC(j,j)
               se(nef+j)=1.d0
            else
            do i=1,j
               ij=i+j*(j-1)/2
               b(nef+ij)=VC(i,j)
               se(nef+ij)=1.d0
            end do
            end if
         end do

      end if

 1589 continue

      deallocate(Y,X,idprob,idea,idg,nmes)

      return
      end subroutine hetmixlin



C-----------------------------------------------------------
C                        FUNCPA
C------------------------------------------------------------


      double precision function funcpa(b,npm,id,thi,jd,thj)



      use commun



      implicit none
      integer ::i,j,k,l,m,g,l2,m2
      integer ::id,jd,jj
      integer ::npm,nef
      integer ::ier,nmoins

      double precision,dimension(maxmes,nv) ::Z,P,X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: vrais,eps,det
      double precision ::thi,thj,temp
      double precision ::Y4,expo

      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3
      double precision,dimension(ng) :: pi

      b1=0.d0
      eps=1.D-20
      do k=1,npm
         b1(k)=b(k)
      end do

      if (id.ne.0) b1(id)=b1(id)+thi
      if (jd.ne.0) b1(jd)=b1(jd)+thj

c----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng

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

c ----------- boucle sur les individus -------------

      vrais=0.d0
      do i=1,ns


c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(i,j,k))
               end do
            end if

         end do



c creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0
         do j=1,nmes(i)
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(i,j))
         end do

c creation de P=Zi*Ut et V=P*P' que si non spec aux classes

         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Se

c Vi en vecteur
            jj=0
            Vi=0.d0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL DSINV(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               funcpa=9999.d0
               goto 654
            end if

c     retransformation du vecteur Vi en matrice :

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
c     debut du calcul de la vraisemblance

       vrais=vrais-nmes(i)*dlog(dble(2*3.14159265))


c contribution individuelle a la vraisemblance

c cas 1 : ng=1

       if (ng.eq.1) then

            vrais=vrais-det
            b0=0.d0
            l=0
            X00=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
                     X00(j,l)=dble(X(i,j,k))
                  end do
                  b0(l)=b1(nprob+l)
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

c cas 2 :  ng>1  composantes
         else

c transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            Xprob=0.d0
            Xprob(1)=1
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(i,1,k)
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

c creation des vecteurs de variables explicatives
            l=0
            m=0
            X00=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
                     X2(j,l)=dble(X(i,j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X00(j,m)=dble(X(i,j,k))
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
                     m2=m2+1
                     b0(m2)=b1(nprob+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nmoins+g)
                     nmoins=nmoins+ng
                  end if
               end do

C variance covariance si spec aux classes :

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

c Vi en vecteur
                  Vi=0.d0
                  jj=0
                  do j=1,nmes(i)
                     do k=j,nmes(i)
                        jj=j+k*(k-1)/2
                        Vi(jj)=VC(j,k)
                     end do
                  end do

                  CALL DSINV(Vi,nmes(i),eps,ier,det)
                  if (ier.eq.-1) then
                     funcpa=9999.d0
                     goto 654
                  end if

c     retransformation du vecteur Vi en matrice :

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
       end do

C FIN BOUCLE SUJET

      funcpa=vrais/2.D0      

 654  continue

      return


      end

C-------------------------------------------------------------
C                   MARQ98
C-------------------------------------------------------------


      subroutine marq98(b,m,ni,v,rl,ier,istop,ca,cb,dd)

c
c  fu = matrice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion matrice d'info (ier=1)

c  5: vraisemblance positive

      use parameters

      IMPLICIT NONE
      integer m,ni,nql,ii,nfmax,idpos
      integer ier,istop,igrad,ncount,id,jd
      integer m1,j,i,ij,k
      integer npl
      double precision fu(m*(m+3)/2),da,dm,ga,tr
      double precision b(m),delta(m),ca,cb,rl
      double precision V(m*(m+3)/2),V1(m*(m+3)/2)
      double precision V2(m,m)
      double precision GHG
      double precision funcpa,det, step,eps
      double precision b1(m),vw,fi,maxt,z,bh(m),rl1,th,ep,dd


C ----------------------------------------------------


      npl=2
      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-6
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
 10   continue

      call deriva(b,m,v,rl)

      rl1=rl

      dd = 0.d0
      do 13 i=m*(m+1)/2+1,m*(m+3)/2
         dd = dd + v(i)*v(i)
 13   continue
      dd=dd/dabs(RL)

      do i=1,m
         do j=i,m
            ij=(j-1)*j/2+i
            fu(ij)=v(ij)
            V2(i,j)=v(ij)
            v2(j,i)=v(ij)
         end do
      END DO
      call dsinv(fu,m,ep,ier,det)
      if (ier.eq.-1) then
         dd=1.d0
      else
         GHG = 0.d0
         do 14 i=1,m
            do 14 j=1,m
               if(j.ge.i) then
                  ij=(j-1)*j/2+i
               else
                  ij=(i-1)*i/2+j
               end if
               GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
 14         continue
         dd=GHG/dble(m)
      end if

      if ( ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100

      tr=0.d0
      do 300 i=1,m
         ii=i*(i+1)/2
         tr=tr+dabs(v(ii))
 300  continue
      tr=tr/dble(m)

      ncount=0
      ga=0.01d0
 400  do 1 i=1,nfmax+m
         fu(i)=v(i)
 1    continue
      do 500 i=1,m
         ii=i*(i+1)/2
         if (v(ii).ne.0) then
            fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
         else
            fu(ii)=da*ga*tr
         endif
 500  continue
      call dchole(fu,m,nql,idpos)

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
         do 2 i=1,m
            delta(i)=fu(nfmax+i)
            b1(i)=b(i)+delta(i)
 2       continue
         rl=funcpa(b1,m,id,z,jd,z)
         if(rl.eq.9999.d0) then
            istop=4
            goto 110
         end if
         igrad=1
         if (rl1.lt.rl) then
            if(da.lt.eps) then
               da=eps
            else
               da=da/(dm+2.d0)
            endif
            goto 800
         endif
      endif
C      write(6,*) 'loglikelihood not improved '
      if(dabs(maxt(delta,m)).eq.0.D0) then
         vw=1.D-5
      else
         vw=th/maxt(delta,m)
      endif
      step=dlog(1.5d0)
C      write(*,*) 'searpas'
      call searpas(vw,step,b,bh,m,delta,fi)
      rl=-fi
      do 4 i=1,m
         delta(i)=vw*delta(i)
 4    continue
      da=(dm-3.d0)*da

 800  cb=dabs(rl1-rl)
      ca=0.d0
      do 5 i=1,m
         ca=ca+delta(i)*delta(i)
 5    continue
C      write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
      do 6 i=1,m
         b(i)=b(i)+delta(i)
 6    continue

      ni=ni+1
      if (ni.ge.70) then
         istop=2
C        write(6,*) 'maximum number of iteration reached'
         goto 110
      end if
      goto 10
c********************
c     inversion matrice d'information
c********************
 100  continue
      istop=1
      DO k=1,npl*(npl+1)/2
         v1(k)=v(k)
      END DO
C      write(*,*)'Final ca,cb,dd',ca,cb,dd
      ep=10.d-10

      call dsinv(v,m,ep,ier,det)
      if (ier.eq.-1) then
         istop=3
         call dsinv(v1,npl,ep,ier,det)
         if (ier.eq.-1) then
            istop=31
         else
            DO k=1,npl*(npl+1)/2
               v(k)=v1(k)
            END DO
         end if
      endif

 110  continue
      return
      end subroutine marq98



C------------------------------------------------------------
C                          DERIVA
C------------------------------------------------------------

      subroutine deriva(b,m,v,rl)


      implicit none
      integer i0,m,m1,ll,i,k,j
      double precision funcpa,thn,th,z,rl,b(m),vl,temp
      double precision fcith(m),th2,thi,thj
      double precision V(m*(m+3)/2)

c
c     v:matrice d'information+score
c     calcul de la derivee premiere par "central difference"
c     calcul des derivees secondes par "forward difference"
c

      th=1.d-5
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      rl=funcpa(b,m, i0,z,i0,z)

      if(rl.eq.9999.d0) then
         go to 123
      end if

C      write(*,*)'rl dans deriva',rl
      do 2 i=1,m

         th=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))

         fcith(i)=funcpa(b,m,i,th,i0,z)
         if(fcith(i).eq.9999.d0) then
            rl=9999.d0
            go to 123
         end if

 2    continue
      k=0
      m1=m*(m+1)/2
      ll=m1
      do i=1,m
         ll=ll+1
         thn=-DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
         temp=funcpa( b,m, i,thn,i0,z)
         if(temp.eq.9999.d0) then
            rl=9999.d0
            goto 123
         end if
         vl=(fcith(i)-temp)/(2.d0*(-thn))
         v(ll)=vl
         do 1 j=1,i
            k=k+1

            thi=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
            thj=DMAX1(1.d-7, 1.d-4 * DABS(b(j)))

            temp=funcpa(b,m,i,thi,j,thj)
            if(temp.eq.9999.d0) then
               rl=9999.d0
               goto 123
            end if
            v(k)=-(temp-fcith(j)-fcith(i)+rl)
     &           /(thi*thj)
 1       continue
      end do
 123   continue
      return
      end subroutine deriva



C------------------------------------------------------------
C                        SEARPAS
C------------------------------------------------------------


      SUBROUTINE SEARPAS(VW,STEP,B,BH,M,DELTA,FIM)
C
C  MINIMISATION UNIDIMENSIONNELLE
C
      implicit none

       INTEGER I,M
       DOUBLE PRECISION VLW,VLW1,VLW2,VLW3,VW,VM,B(m)
       DOUBLE PRECISION FI1,FI2,FI3,FIM,BH(m),DELTA(m)
       DOUBLE PRECISION STEP
C
C
       VLW1=DLOG(VW)
       VLW2=VLW1+STEP
       CALL VALFPA(VLW1,FI1,B,BH,M,DELTA)
       CALL VALFPA(VLW2,FI2,B,BH,M,DELTA)
C
       IF(FI2.GE.FI1) THEN
          VLW3=VLW2
          VLW2=VLW1
          FI3=FI2
          FI2=FI1
C
          STEP=-STEP
C
          VLW1=VLW2+STEP
          CALL VALFPA(VLW1,FI1,B,BH,M,DELTA)

          IF (FI1.GT.FI2) GOTO 50
       ELSE
          VLW=VLW1
          VLW1=VLW2
          VLW2=VLW
          FIM=FI1
          FI1=FI2
          FI2=FIM
       ENDIF
C
       DO 20 I=1,40
          VLW3=VLW2
          VLW2=VLW1
          FI3=FI2
          FI2=FI1
C
          VLW1=VLW2+STEP
          CALL VALFPA(VLW1,FI1,B,BH,M,DELTA)
          IF(FI1.GT.FI2) GO TO 50
          IF (FI1.eq.FI2) THEN
             FIM=FI2
             VM=VLW2
             GO TO 100
          ENDIF
 20    CONTINUE
C
C  PHASE 2 APPROXIMATION PAR QUADRIQUE
C
50     CONTINUE
C
C  CALCUL MINIMUM QUADRIQUE
C
         VM=VLW2-STEP*(FI1-FI3)/(2.d0*(FI1-2.d0*FI2+FI3))
         CALL VALFPA(VM,FIM,B,BH,M, DELTA)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
C      write(*,*)'fi1,fi2,fi3,fim',fi1,fi2,fi3,fim
      VW=DEXP(VM)
      RETURN

      END subroutine searpas


C------------------------------------------------------------
C                          VALFPA
C------------------------------------------------------------

        subroutine valfpa(vw,fi,b,bk,m,delta)

        implicit none
        integer m,i,i0
        double precision vw,fi,delta(m),bk(m)
        double precision funcpa,b(m),z
         z=0.d0
         i0=1
         do 1 i=1,m
            bk(i)=b(i)+dexp(vw)*delta(i)
1        continue
         fi=-funcpa(bk,m,i0,z,i0,z)

         return
         end subroutine valfpa

C------------------------------------------------------------
C                            MAXT
C------------------------------------------------------------


      DOUBLE PRECISION FUNCTION MAXT(DELTA,M)


      implicit none


      INTEGER I,M
      DOUBLE PRECISION DELTA(m)
c
      MAXT=Dabs(DELTA(1))
      DO 2 I=2,M
      IF(DABS(DELTA(I)).GT.MAXT) MAXT=DABS(DELTA(I))
2     CONTINUE
      RETURN
      END


C------------------------------------------------------------
C                         DCHOLE
C------------------------------------------------------------

      subroutine dchole(a,k,nq,idpos)
      implicit none
      double precision a(1000),term,xn,diag,p
      dimension is(200)
      equivalence (term,xn)
      integer k,nq,i,ii,i1,i2,i3,m,is,j,k2,jmk
      integer ijm,irm,jji,jjj,l,jj,iil,jjl,il,idpos
c
c
c      ss programme de resolution d'un systeme lineaire symetrique
c
c
c       k ordre du systeme /
c       nq nombre de seconds membres
c
c       en sortie les seconds membres sont remplaces par les solutions
c       correspondantes
c

      i2=0
      ii=0
      idpos=0
      k2=k+nq
c     calcul des elements de la matrice
      do 13 i=1,k
      ii=i*(i+1)/2
c       elements diagonaux
      diag=a(ii)
      i1=ii-i
      if(i-1.ne.0) goto 1
      if(i-1.eq.0) goto 4
1     i2=i-1
      do 3 l=1,i2
      m=i1+l
      p=a(m)
      p=p*p
      if(is(l).lt.0) goto 2
      if(is(l).ge.0) goto 3
2     p=-p
3     diag=diag-p
4     if(diag.lt.0) goto 5
      if(diag.eq.0) goto 50
      if(diag.gt.0) goto 6
5     is(i)=-1
      idpos=idpos+1
      diag=-dsqrt(-diag)
      a(ii)=-diag
      go to 7
6     is(i)=1
      diag=dsqrt(diag)
      a(ii)=diag
c       elements non diagonaux
7     i3=i+1
      do 13 j=i3,k2
      jj=j*(j-1)/2+i
      jmk=j-k-1
      if(jmk.le.0) goto 9
      if(jmk.gt.0) goto 8
8     jj=jj-jmk*(jmk+1)/2
9     term=a(jj)
      if(i-1.ne.0) goto 10
      if(i-1.eq.0) goto 13
10    do 12 l=1,i2
      iil=ii-l
      jjl=jj-l
      p=a(iil)*a(jjl)
      il=i-l
      if(is(il).lt.0) goto 11
      if(is(il).ge.0) goto 12
11    p=-p
12    term=term-p
13    a(jj)=term/diag
c       calcul des solutions
      jj=ii-k+1
      do 45 l=1,nq
      jj=jj+k
      i=k-1
14    jji=jj+i
      xn=a(jji)
      if(i-k+1.lt.0) goto 20
      if(i-k+1.ge.0) goto 22
20    j=k-1
21    jjj=jj+j
      ijm=i+1+j*(j+1)/2
      xn=xn-a(jjj)*a(ijm)
      if(j-i-1.le.0) goto 22
      if(j-i-1.gt.0) goto 30
30    j=j-1
      go to 21
22    irm=(i+1)*(i+2)/2
      a(jji)=xn/a(irm)
      if(i.le.0) goto 45
      if(i.gt.0) goto 40
40    i=i-1
      go to 14
45    continue
50    continue
      return
      end subroutine dchole


C------------------------------------------------------------
C                               DMFSD
C------------------------------------------------------------


      SUBROUTINE DMFSD(A,N,EPS,IER)
C
C   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
C   MATRICE = TRANSPOSEE(T)*T
C   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
C            PAR COLONNE DE LA METRICE A FACTORISER
C   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
C
C   SUBROUTINE APPELE PAR DSINV
C
C   N : DIM. MATRICE
C   EPS : SEUIL DE TOLERANCE
C   IER = 0 PAS D'ERREUR
C   IER = -1 ERREUR
C   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
C
      implicit none
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION DPIV,DSUM
      DOUBLE PRECISION EPS,TOL
      INTEGER I,K,L,N,IER,KPIV,IND,LEND,LANF,LIND
C
C   TEST ON WRONG INPUT PARAMETER N
C
      DPIV=0.d0
      IF(N-1.lt.0) goto 12
      IF(N-1.ge.0) goto 1
1     IER=0
C
C   INITIALIZE DIAGONAL-LOOP
C
      KPIV=0
      DO 11 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
C
C   CALCULATE TOLERANCE
C
      TOL=DABS(EPS*SNGL(A(KPIV)))
C
C   START FACTORIZATION-LOOP OVER K-TH ROW
C
      DO 11 I=K,N
      DSUM=0.D0
      IF(LEND.lt.0) goto 2
      IF(LEND.eq.0) goto 4
      IF(LEND.gt.0) goto 2
C
C   START INNER LOOP
C
2     DO 3 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
3     DSUM=DSUM+A(LANF)*A(LIND)
C
C   END OF INNEF LOOP
C
C   TRANSFORM ELEMENT A(IND)
C
4     DSUM=A(IND)-DSUM
      IF(I-K.ne.0) goto 10
      IF(I-K.eq.0) goto 5
C
C   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
C
5     IF(SNGL(DSUM)-TOL.le.0) goto 6
      IF(SNGL(DSUM)-TOL.gt.0) goto 9
6     IF(DSUM.le.0) goto 12
      IF(DSUM.gt.0) goto 7
7     IF(IER.le.0) goto 8
      IF(IER.gt.0) goto 9
8     IER=K-1
C
C   COMPUTE PIVOT ELEMENT
C
9     DPIV=DSQRT(DSUM)
      A(KPIV)=DPIV
      DPIV=1.D0/DPIV
      GO TO 11
C
C   CALCULATE TERMS IN ROW
C
10    A(IND)=DSUM*DPIV
11    IND=IND+I
C
C   END OF DIAGONAL-LOOP
C
      RETURN
12    IER=-1
      RETURN

      END subroutine dmfsd


C------------------------------------------------------------
C                            DSINV
C------------------------------------------------------------


      SUBROUTINE DSINV(A,N,EPS,IER,DET)
C
C     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
C
C     MATRICE = TRANSPOSEE(T)*T
C     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
C
C     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
C         STOCKEE COLONNE PAR COLONNE
C     DIM. MATRICE A INVERSER = N
C     DIM. TABLEAU A = N*(N+1)/2
C
C     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
C           COMME NUL
C
C     IER : CODE D'ERREUR
C         IER=0 PAS D'ERREUR
C         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
C         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
C
      implicit none
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION DIN,WORK,DET
      DOUBLE PRECISION EPS
      INTEGER N,IER,IND,IPIV,I,J,K,L,MIN,KEND
      INTEGER LHOR,LVER,LANF
C
C     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
C     A=TRANSPOSE(T) * T
C
      CALL DMFSD(A,N,EPS,IER)
      IF(IER.lt.0) goto 9
      IF(IER.ge.0) goto 1
C
C     INVERT UPPER TRIANGULAR MATRIX T
C     PREPARE INVERSION-LOOP
C
C
CC calcul du log du determinant
1     DET=0.D0
      DO 10 i=1,N
        DET=DET+dlog(A(i*(i+1)/2))
10    CONTINUE
      DET=2*DET

      IPIV=N*(N+1)/2
      IND=IPIV
C
C     INITIALIZE INVERSION-LOOP
C
      DO 6 I=1,N
      DIN=1.D0/A(IPIV)
      A(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF(KEND.le.0) goto 5
      IF(KEND.gt.0) goto 2
2     J=IND
C
C     INITIALIZE ROW-LOOP
C
      DO 4 K=1,KEND
      WORK=0.D0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
C
C     START INNER LOOP
C
      DO 3 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
3     WORK=WORK+A(LVER)*A(LHOR)
C
C     END OF INNER LOOP
C
      A(J)=-WORK*DIN
4     J=J-MIN
C
C     END OF ROW-LOOP
C
5     IPIV=IPIV-MIN
6     IND=IND-1
C
C     END OF INVERSION-LOOP
C
C     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
C     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
C     INITIALIZE MULTIPLICATION-LOOP
C
      DO 8 I=1,N
      IPIV=IPIV+I
      J=IPIV
C
C     INITIALIZE ROW-LOOP
C
      DO 8 K=I,N
      WORK=0.D0
      LHOR=J
C
C     START INNER LOOP
C
      DO 7 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+A(LHOR)*A(LVER)
7     LHOR=LHOR+L
C
C     END OF INNER LOOP
C
      A(J)=WORK
8     J=J+K
C
C     END OF ROW-AND MULTIPLICATION-LOOP
C
9     RETURN

      END SUBROUTINE dsinv


C------------------------------------------------------------
C                      POSTPROB
C------------------------------------------------------------

c-------------------------------------------------------------
c
c          Subroutine pour calculer les
c      probabilites a posteriori de suivre chacune
c      des composantes g pour chacun des sujets i
c
c                                 25/02/03
c-------------------------------------------------------------

      subroutine postprob(b,npm,PPI)


      use commun


      implicit none

	  integer ::i,j,k,l,m,g,l2,m2
      integer ::jj
      integer ::npm,nef
      integer ::ier,nmoins

      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision ::Y4,f
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3



      eps=1.D-20

      PPI=0.D0

      do k=1,npm
         b1(k)=b(k)
      end do


c----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng
      
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
      

c ----------- boucle sur les individus -------------

      do i=1,ns

c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(i,j,k))
               end do
            end if

         end do



c creation de s2*I et Y1

         Se=0.d0
		 Y1=0.d0
         do j=1,nmes(i)
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(i,j))
         end do

c creation de P=Zi*Ut et V=P*P' que si non spec aux classes

        if (nwg.eq.0) then
          P=0.d0
          P=MATMUL(Z,Ut)
          VC=0.d0
          VC=MATMUL(P,transpose(P))+Se

c Vi en vecteur
         jj=0
         do j=1,nmes(i)
            do k=j,nmes(i)
               jj=j+k*(k-1)/2
               Vi(jj)=VC(j,k)
            end do
         end do

         CALL DSINV(Vi,nmes(i),eps,ier,det)
         if (ier.eq.-1) then
            PPI=-1.d0
            go to 147
         end if

c retransformation du vecteur Vi en matrice :

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


c transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       Xprob=0.d0
       Xprob(1)=1
       l=0
       do k=1,nv
          if (idprob(k).eq.1) then
             l=l+1
             Xprob(1+l)=X(i,1,k)
          end if
       end do
C     write(*,*)'l apres Xprob',l,(Xprob(j),j=1,10)
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
       
C     write(*,*)'pi',(pi(g),g=1,ng)
       
c     creation des vecteurs de variables explicatives
       l=0
       m=0
       X0=0.d0
       X2=0.d0
       do k=1,nv
          if (idg(k).eq.2) then
             l=l+1
             do j=1,nmes(i)
                X2(j,l)=dble(X(i,j,k))
             end do
          else if (idg(k).eq.1) then
             m=m+1
             do j=1,nmes(i)
                X0(j,m)=dble(X(i,j,k))
             end do
          end if
       end do
       
       
c     calcul de la vraisemblance par composante
       f=0.d0
       fi=0.d0
       do g=1,ng
          nmoins=0
          l2=0
          m2=0
          do k=1,nv
             if (idg(k).eq.1) then
                m2=m2+1
                b0(m2)=b1(nprob+nmoins+1)
                nmoins=nmoins+1
             else if (idg(k).eq.2) then
                l2=l2+1
                b2(l2)=b1(nprob+nmoins+g)
                nmoins=nmoins+ng
             end if
          end do
          
C     variance covariance si spec aux classes :
          
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
             
c     Vi en vecteur
             jj=0
             do j=1,nmes(i)
                do k=j,nmes(i)
                   jj=j+k*(k-1)/2
                   Vi(jj)=VC(j,k)
                end do
             end do
             
             CALL DSINV(Vi,nmes(i),eps,ier,det)
             if (ier.eq.-1) then
                PPI=-1.d0
                goto 147
             end if
             
c     retransformation du vecteur Vi en matrice :
             
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
       
      end do
      
 147  continue
      return
      
      end subroutine postprob
      
C------------------------------------------------------------
C                      RESIDUALS
C------------------------------------------------------------


      subroutine residuals(b1,npm,ppi,resid_m,pred_m_g,resid_ss,
     &     pred_ss_g)

      use commun

      implicit none
      integer ::i,j,k,l,m,g,l2,m2
      integer ::jj
      integer ::npm,nef
      integer ::ier,nmoins,nmes_cur,n2,nmoins2

      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob,err2
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se,VC1
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob,b3
      double precision :: eps,det
      double precision ::temp
      double precision,dimension(nv,maxmes)::Valea

      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m
     &     ,pred_m,resid_ss,pred_ss
      double precision,dimension(nobs,ng)::
     &     pred_m_g,pred_ss_g

      double precision,dimension(ns,ng) ::PPI






      eps=1.D-20

c----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng



C      write(*,*)'nvc',nvc,'nea',nea,'nwg',nwg,'nef',nef
c creation de Ut, decomposition de cholesky pour G
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

C         write(*,*)'Ut',Ut

c ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m=0.d0
      pred_m_g=0.d0
      pred_ss=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0

      nmes_cur=0
      do i=1,ns

c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(i,j,k))
               end do
            end if

         end do



c     creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0
         do j=1,nmes(i)
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(i,j))
         end do

c     creation de P=Zi*Ut et V=P*P' que si non spec aux classes


c     debut du calcul de la vraisemblance

c     cas 1 : ng=1

         if (ng.eq.1) then


            Valea=0.d0
            VC=0.d0
            P=0.d0


            P=MATMUL(Z,Ut)
            Valea=MATMUL(Ut,transpose(P))
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Se


c     Vi en vecteur

            jj=0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL DSINV(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               do j=1,nmes(i)
                  resid_m(nmes_cur+j)=9999.d0
                  pred_m(nmes_cur+j)=9999.d0
                  pred_m_g(nmes_cur+j,1)=9999.d0

                  resid_ss(nmes_cur+j)=9999.d0
                  pred_ss(nmes_cur+j)=9999.d0
                  pred_ss_g(nmes_cur+j,1)=9999.d0
               end do
               goto 654
            end if

c     retransformation du vecteur Vi en matrice :
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
                     X0(j,l)=dble(X(i,j,k))
                  end do
                  b0(l)=b1(nprob+l)
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
               pred_m(nmes_cur+j)=mu(j)
               pred_m_g(nmes_cur+j,1)=mu(j)

               resid_ss(nmes_cur+j)=Y2(j)-err1(j)
               pred_ss(nmes_cur+j)=pred1(j)
               pred_ss_g(nmes_cur+j,1)=pred_ss(nmes_cur+j)

            end do


c     cas 2 :  ng>1  composantes
         else


c     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            Xprob=0.d0
            Xprob(1)=1.d0
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(i,1,k)
               end if
            end do
C     write(*,*)'l apres Xprob',l,(Xprob(j),j=1,10)
            pi=0.d0
            temp=0.d0
            Do g=1,ng-1
               bprob=0.d0
               do k=1,nvarprob
                  bprob(k)=b1((k-1)*(ng-1)+g)
               end do

C     write(*,*)'g=',g,'nvarprob',nvarprob,'bprob='
C     &,(bprob(k),k=1,nvarprob)

               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

            end do
            pi(ng)=1.d0/(1.d0+temp)
            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do

C     write(*,*)'pi',(pi(g),g=1,ng)

c     creation des vecteurs de variables explicatives
            l=0
            m=0
            X0=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
                     X2(j,l)=dble(X(i,j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X0(j,m)=dble(X(i,j,k))
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
                     m2=m2+1
                     b0(m2)=b1(nprob+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nmoins+g)
                     nmoins=nmoins+ng
                  end if
                  IF (IDEA(k).EQ.1.and.idg(k).eq.1) THEN
                     n2=n2+1
                     b3(n2)=b1(nprob+nmoins2+1)
                     nmoins2=nmoins2+1
                  else if(IDEA(k).EQ.1.and.idg(k).eq.2) THEN
                     n2=n2+1
                     b3(n2)=b1(nprob+nmoins2+g)
                     nmoins2=nmoins2+1
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


c     Vi en vecteur

               jj=0
               do j=1,nmes(i)
                  do k=j,nmes(i)
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do
               CALL DSINV(Vi,nmes(i),eps,ier,det)
               if (ier.eq.-1) then
                  do j=1,nmes(i)
                     resid_m(nmes_cur+j)=9999.d0
                     pred_m(nmes_cur+j)=9999.d0

                     resid_ss(nmes_cur+j)=9999.d0
                     pred_ss(nmes_cur+j)=9999.d0
                     do l=1,ng
                        pred_ss_g(nmes_cur+j,l)=9999.d0
                        pred_m_g(nmes_cur+j,l)=9999.d0
                     end do
                  end do
                  goto 654
               end if

c     retransformation du vecteur Vi en matrice :
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
                  pred_m_g(nmes_cur+j,g)=mu(j)
                  pred_ss_g(nmes_cur+j,g)=pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j)
     &                 +ppi(i,g)*(Y1(j)-pred1(j))
                  pred_ss(nmes_cur+j)=pred_ss(nmes_cur+j)+ppi(i,g)
     &                 *pred1(j)
                  pred_m(nmes_cur+j)=pred_m(nmes_cur+j)+pi(g)*mu(j)
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))

               end do
            end do

         end if


         nmes_cur=nmes_cur+nmes(i)
      end do

C FIN BOUCLE SUJET

 654  continue

      end subroutine residuals

