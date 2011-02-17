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
c       C cile Proust, INSERM U897, ISPED,
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
c maxiter : number of iteration
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

      module commun

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg
     & ,nprob,nvarprob,maxmes,nobs
      double precision,dimension(:),allocatable,save::Y
      double precision,dimension(:,:),allocatable,save ::X
      integer,dimension(:),allocatable,save ::idea,idg,idprob
      integer,dimension(:),allocatable,save :: nmes,prior
      end module commun

c*************************************************************
      module parameters
          double precision,save::epsa,epsb,epsd
          integer,save::maxiter
      end module parameters




C================ SUBROUTINES ================================

      subroutine hetmixlin(
     &      Y0,X0,Prior0,
     &      idprob0,idea0,idg0
     &     ,ns0,ng0,nv0,nobs0
     &     ,nea0
     &     ,nmes0,idiag0,nwg0
     &     ,npm0,b,vopt
     &     ,vrais
     &     ,ni,istop,gconv
     &     ,ppi0
     &     ,resid_m0,resid_ss0
     &     ,pred_m_g0
     &     ,pred_ss_g0
     &     ,pred_RE
     &     ,convB,convL,convG,maxiter0)

      use commun
      use parameters
C      use optim

      IMPLICIT NONE


	!D claration des variables en entree
      integer,intent(in):: nv0,maxiter0,nea0
      integer, intent(in) :: ns0, ng0, nobs0, idiag0, nwg0, npm0
      integer, dimension(nv0), intent(in) :: idea0,idg0,idprob0
      integer, dimension(ns0), intent(in) :: nmes0,Prior0
      double precision, dimension(nobs0), intent(in) :: Y0
      double precision, dimension(nobs0*nv0), intent(in) :: X0
      double precision, intent(in) :: convB, convL, convG
	!D claration des variable en entr e et sortie
      double precision, dimension(npm0), intent(inout) :: b

	!D claration des variables en sortie
      double precision, intent(out) :: vrais
      double precision, dimension(3), intent(out) :: gconv
      double precision, dimension(ns0*ng0), intent(out) :: ppi0
      double precision, dimension(nobs0), intent(out) :: resid_m0
      double precision, dimension(nobs0), intent(out) :: resid_ss0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_m_g0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_ss_g0
      double precision, dimension(npm0*(npm0+1)/2),intent(out) :: vopt
      double precision, dimension(npm0*(npm0+3)/2) :: v
      integer, intent(out) :: ni, istop
	!Variables locales
      integer :: jtemp,nef,i,g,j,ij,npm,ier,k,ktemp,ig,nmestot,it
      double precision :: eps, ca, cb, dd
      double precision, dimension(ns0,ng0) :: PPI
      double precision, dimension(npm0) :: mvc
      double precision, dimension(ns0*nea0), intent(out)::pred_RE
      double precision, dimension(nobs0) :: resid_m, resid_ss
      double precision, dimension(nobs0,ng0):: pred_m_g, pred_ss_g


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

      allocate(Y(ns0*maxmes),idprob(nv0),X(ns0*maxmes,nv0)
     & ,idea(nv0),idg(nv0),nmes(ns0),prior(ns0))



      ppi0=1.d0
      resid_m0=0.d0
      resid_ss0=0.d0
      pred_m_g0=0.d0
      pred_ss_g0=0.d0
      pred_re=0.d0

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
      prior=0
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

C creation des parametres

      nea=nea0
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
      if (nwg.gt.0) then
         do i=1,nwg
            B(nef+nvc+i)=abs(B(nef+nvc+i))
         end do
      end if

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

         do k=1,nwg+1
            b(nef+nvc+k)=abs(b(nef+nvc+k))
         end do

         vopt=0.d0
         vopt(1:npm0*(npm0+1)/2)=v(1:npm0*(npm0+1)/2)

C probas posteriori


         if (istop.eq.1) then
            if (ng.gt.1) then
               call postprob(B,npm,PPI)
            end if

            call residuals(b,npm,ppi,resid_m,pred_m_g,resid_ss
     &           ,pred_ss_g,pred_RE)

            ig=0
            ij=0
            do i=1,ns
                if (ng.gt.1) then
                    do g=1,ng0
                        ig=ig+1
                        ppi0(ig)=PPI(i,g)
                    end do
                end if
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

      end if

 1589 continue

      deallocate(Y,X,idprob,idea,idg,nmes,prior)

      return
      end subroutine hetmixlin



C-----------------------------------------------------------
C                        FUNCPA
C------------------------------------------------------------


      double precision function funcpa(b,npm,id,thi,jd,thj)

      use commun
C      use optim

      IMPLICIT NONE

      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,nef,it
      integer ::ier,nmoins,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
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
       kk=0
       vrais=0.d0
       it=0
       do i=1,ns
         if (i==1) then
            it=0
         else
            it=it+nmes(i-1)
         end if

c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

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



c creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0
         do j=1,nmes(i)
            kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(kk))
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
               funcpa=-1.d9
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
                        X00(j,l)=dble(X(it+j,k))
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
            if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
            else
c transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
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
c creation des vecteurs de variables explicatives
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
                     funcpa=-1.d9
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

            if (expo.le.0.d0) then
                funcpa=-1.d9
                goto 654
            end if
            vrais=vrais+2*log(expo)
         end if
      end do

C FIN BOUCLE SUJET

      funcpa=vrais/2.D0

 654  continue

      return


      end function funcpa


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
C      use optim


      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,nef,ier,nmoins,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
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
      kk=0
	it=0
      do i=1,ns
	if (i==1) then
		it=0
	else
		it=it+nmes(i-1)
	end if

c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

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



c creation de s2*I et Y1

         Se=0.d0
		 Y1=0.d0
         do j=1,nmes(i)
	     kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
	     Y1(j)=dble(Y(kk))
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

       if (prior(i).ne.0) then
          pi=0.d0
          pi(prior(i))=1.d0
       else
c transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       Xprob=0.d0
       Xprob(1)=1
       l=0
       do k=1,nv
          if (idprob(k).eq.1) then
             l=l+1
             Xprob(1+l)=X(it+1,k)
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
      endif
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
		X2(j,l)=dble(X(it+j,k))
             end do
          else if (idg(k).eq.1) then
             m=m+1
             do j=1,nmes(i)
                X0(j,m)=dble(X(it+j,k))
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
     &     pred_ss_g,pred_RE)

      use commun
C      use optim

      implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,nef
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea) :: err2
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se,VC1
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob,b3
      double precision :: eps,det
      double precision ::temp
      double precision,dimension(nea,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m
     &     ,pred_m,resid_ss,pred_ss
      double precision,dimension(nobs,ng)::
     &     pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(ns*nea)::pred_RE






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
      pred_re=0.d0
      nmes_cur=0
      kk=0
      it=0

       do i=1,ns

	if (i==1) then
		it=0
	else
		it=it+nmes(i)
	end if

c -------- creation de Vi = ZiGZi'+se*seIni ----------

c creation de Zi

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



c     creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0
         do j=1,nmes(i)
            kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(kk))
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
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
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
                     X0(j,l)=dble(X(it+j,k))
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


            do k=1,nea
               pred_RE((i-1)*nea+k)=err2(k)
            end do


c     cas 2 :  ng>1  composantes
         else

            if (prior(i).ne.0) then
	           pi=0.d0
               pi(prior(i))=1.d0
            else
c     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            Xprob=0.d0
            Xprob(1)=1.d0
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(it+1,k)
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
          end if
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
                     X2(j,l)=dble(X(it+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X0(j,m)=dble(X(it+j,k))
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
                  do k=1,nea
                    pred_RE((i-1)*nea+k)=99999.d0
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

               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*
     &                 err2(k)
               end do

            end do

         end if

 654  continue

          nmes_cur=nmes_cur+nmes(i)
       enddo

C FIN BOUCLE SUJET


      end subroutine residuals







C ================== MODULE OPTIM2011=========================



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
c  4: Erreur

      use parameters

      IMPLICIT NONE

C   variables globales
      integer,intent(in) :: m
      integer,intent(inout)::ni,ier,istop
      double precision,dimension(m*(m+3)/2),intent(out)::v
      double precision,intent(out)::rl
      double precision,dimension(m),intent(inout)::b
      double precision,intent(inout)::ca,cb,dd

C   variables locales
      integer::nql,ii,nfmax,idpos,ncount,id,jd
     & ,m1,j,i,ij
      double precision,dimension(m*(m+3)/2)::fu
      double precision,dimension(m)::delta,b1,bh
      double precision::da,dm,ga,tr
      double precision::GHG,funcpa,det,step,eps
     & ,vw,fi,z,rl1,th,ep
      double precision::maxt

C ----------------------------------------------------

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

      Main:Do
        call deriva(b,m,v,rl)
        rl1=rl

        dd = 0.d0
        do i=m*(m+1)/2+1,m*(m+3)/2
           dd = dd + v(i)*v(i)
        end do
        dd=dd/dabs(RL)
c---------
c---------
        fu=0.D0
        do i=1,m
           do j=i,m
              ij=(j-1)*j/2+i
              fu(ij)=v(ij)
           end do
        end do
c---------
c---------- cet appel ne passe pas ----------------- !

        call dsinv(fu,m,ep,ier,det)

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
            do i=1,m
               delta(i)=fu(nfmax+i)
               b1(i)=b(i)+delta(i)
            end do
            rl=funcpa(b1,m,id,z,jd,z)
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
         if(maxt(delta,m).eq.0.D0) then
            vw=th
         else
            vw=th/maxt(delta,m)
         endif
         step=dlog(1.5d0)
C      write(*,*) 'searpas'
         call searpas(vw,step,b,bh,m,delta,fi)
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
C        write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
         do i=1,m
            b(i)=b(i)+delta(i)
         end do

         ni=ni+1
         if (ni.ge.maxiter) then
            istop=2
            write(6,*) 'maximum number of iteration reached'
            goto 110
         end if
      End do Main

      istop=1

      v(1:(m*(m+1)/2))=fu(1:(m*(m+1)/2))



 110   continue
       return
       end subroutine marq98

C------------------------------------------------------------
C                          DERIVA
C------------------------------------------------------------

      subroutine deriva(b,m,v,rl)

      implicit none

      integer,intent(in)::m
      double precision,intent(inout)::rl
      double precision,dimension(m),intent(in)::b
      double precision,dimension((m*(m+3)/2)),intent(out)::v
      double precision,dimension(m)::fcith
      integer ::i0,m1,ll,i,k,j
      double precision::funcpa,thn,th,z,vl,temp,thi,thj
c
c     v:matrice d'information+score
c     calcul de la derivee premiere par "central difference"
c     calcul des derivees secondes par "forward difference"
c

      z=0.d0
      i0=0

      rl=funcpa(b,m,i0,z,i0,z)

      if(rl.eq.-1.d9) then
         goto 123
      end if

      do i=1,m
         th=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
         fcith(i)=funcpa(b,m,i,th,i0,z)
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
         temp=funcpa( b,m, i,thn,i0,z)
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

            temp=funcpa(b,m,i,thi,j,thj)
            if(temp.eq.-1.d9) then
               rl=-1.d9
               exit Main
            end if
            v(k)=-(temp-fcith(j)-fcith(i)+rl)
     &           /(thi*thj)
         end do
      end do Main
 123   continue

      return
      end subroutine deriva
C------------------------------------------------------------
C                        SEARPAS
C------------------------------------------------------------


      subroutine searpas(vw,step,b,bh,m,delta,fim)
C
C  MINIMISATION UNIDIMENSIONNELLE
C
      implicit none

      integer,intent(in)::m
      double precision,dimension(m),intent(in)::b
      double precision,intent(inout)::vw
      double precision,dimension(m),intent(inout)::bh,delta
      double precision,intent(inout)::fim,step
      double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3
      integer::i
C
C
       vlw1=dlog(vw)
       vlw2=vlw1+step
       call valfpa(vlw1,fi1,b,bh,m,delta)
       call valfpa(vlw2,fi2,b,bh,m,delta)
C
       if(fi2.ge.fi1) then
	  vlw3=vlw2
	  vlw2=vlw1
	  fi3=fi2
	  fi2=fi1
	  step=-step
C
          vlw1=vlw2+step
          call valfpa(vlw1,fi1,b,bh,m,delta)
          if(fi1.gt.fi2) goto 50
       else
          vlw=vlw1
          vlw1=vlw2
          vlw2=vlw
          fim=fi1
          fi1=fi2
          fi2=fim
       end if
C
       do i=1,40
          vlw3=vlw2
          vlw2=vlw1
          fi3=fi2
          fi2=fi1
C
          vlw1=vlw2+step
          call valfpa(vlw1,fi1,b,bh,m,delta)
          if(fi1.gt.fi2) goto 50
          if(fi1.eq.fi2) then
             fim=fi2
             vm=vlw2
             goto 100
          end if
       end do
C
C  PHASE 2 APPROXIMATION PAR QUADRIQUE
C
50     continue
C
C  CALCUL MINIMUM QUADRIQUE
C
      vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))
      call valfpa(vm,fim,b,bh,m,delta)
      if(fim.le.fi2) goto 100
      vm=vlw2
      fim=fi2
100   continue
      vw=dexp(vm)

      return

      end subroutine searpas

C------------------------------------------------------------
C                         DCHOLE
C------------------------------------------------------------

      subroutine dchole(a,k,nq,idpos)

      implicit none

      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a

      integer::i,ii,i1,i2,i3,m,j,k2,jmk
      integer::ijm,irm,jji,jjj,l,jj,iil,jjl,il
      integer,dimension(k)::is
      double precision ::term,xn,diag,p
      equivalence (term,xn)


c      ss programme de resolution d'un systeme lineaire symetrique
c
c
c       k ordre du systeme /
c       nq nombre de seconds membres
c
c       en sortie les seconds membres sont remplaces par les solutions
c       correspondantes
c
      jj=0
      i2=0
      ii=0
      idpos=0
      k2=k+nq
c     calcul des elements de la matrice
      do i=1,k
         ii=i*(i+1)/2
c       elements diagonaux
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
c       elements non diagonaux
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
         a(jj)=term/diag
      end do

c       calcul des solutions
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
      end subroutine dchole


      subroutine dmfsd(a,n,eps,ier)
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

      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(inout)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

C
C   TEST ON WRONG INPUT PARAMETER N
C
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
C
C   INITIALIZE DIAGONAL-LOOP
C
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
C
C   CALCULATE TOLERANCE
C
          tol=dabs(eps*sngl(A(kpiv)))
C
C   START FACTORIZATION-LOOP OVER K-TH ROW
C
         do i=k,n
	    dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
C
C   START INNER LOOP
C
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
	       dsum=dsum+A(lanf)*A(lind)
            end do

C
C   END OF INNEF LOOP
C
C   TRANSFORM ELEMENT A(IND)
C
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
C
C   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
C
5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12  !c ici que ça plante  exit do
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
C
C   COMPUTE PIVOT ELEMENT
C
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
C
C   CALCULATE TERMS IN ROW
C
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
	 ind=ind+i
      end do

C
C   END OF DIAGONAL-LOOP
C
      return
12    ier=-1
      return

      end subroutine dmfsd


C------------------------------------------------------------
C                            DSINV
C------------------------------------------------------------


      subroutine dsinv(A,N,EPS,IER,DET)

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

      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(inout)::eps
      double precision,intent(inout),optional::det
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision::din,work
      integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf

C
C     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
C     A=TRANSPOSE(T) * T
C

      call dmfsd(A,n,eps,ier)



      if (ier.lt.0) goto 9
      if (ier.ge.0) det=0.d0
C
C     INVERT UPPER TRIANGULAR MATRIX T
C     PREPARE INVERSION-LOOP
C
C
CC calcul du log du determinant

      do i=1,n
         det=det+dlog(A(i*(i+1)/2))
      end do
      det=2*det
      ipiv=n*(n+1)/2
      ind=ipiv
C
C     INITIALIZE INVERSION-LOOP
C
      do i=1,n
         din=1.d0/A(ipiv)
         A(ipiv)=din
         min=n
         kend=i-1
         lanf=n-kend
         if (kend.le.0) goto 5
         if (kend.gt.0) j=ind
C
C     INITIALIZE ROW-LOOP
C
         do k=1,kend
	    work=0.d0
	    min=min-1
	    lhor=ipiv
	    lver=j
C
C     START INNER LOOP
C
            do l=lanf,min
	        lver=lver+1
		lhor=lhor+l
                work=work+A(lver)*A(lhor)
	    end do
C
C     END OF INNER LOOP
C
            A(j)=-work*din
            j=j-min
	 end do

C
C     END OF ROW-LOOP
C
5        ipiv=ipiv-min
         ind=ind-1
      end do

C
C     END OF INVERSION-LOOP
C
C     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
C     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
C     INITIALIZE MULTIPLICATION-LOOP
C
      do i=1,n
         ipiv=ipiv+i
	 j=ipiv
C
C     INITIALIZE ROW-LOOP
C
	 do k=i,n
	    work=0.d0
	    lhor=j
C
C     START INNER LOOP
C
            do l=k,n
	        lver=lhor+k-i
		work=work+A(lhor)*A(lver)
   		lhor=lhor+l
            end do
C
C     END OF INNER LOOP
C
            A(j)=work
            j=j+k
	 end do
         j=j+k
      end do

C
C     END OF ROW-AND MULTIPLICATION-LOOP
C
9     return
      end subroutine dsinv

C------------------------------------------------------------
C                          VALFPA
C------------------------------------------------------------

        subroutine valfpa(vw,fi,b,bk,m,delta)

        implicit none

        integer,intent(in)::m
        double precision,dimension(m),intent(in)::b,delta
        double precision,dimension(m),intent(out)::bk
        double precision,intent(out)::fi
	double precision::vw,funcpa,z
	integer::i0,i

         z=0.d0
         i0=1
         do i=1,m
            bk(i)=b(i)+dexp(vw)*delta(i)
	 end do
         fi=-funcpa(bk,m,i0,z,i0,z)

         return

         end subroutine valfpa

C------------------------------------------------------------
C                            MAXT
C------------------------------------------------------------


      double precision function maxt(delta,m)

      implicit none

       integer,intent(in)::m
       double precision,dimension(m),intent(in)::delta
       integer::i
c
       maxt=Dabs(delta(1))
       do i=2,m
         if(Dabs(delta(i)).gt.maxt)then
	    maxt=Dabs(delta(i))
	 end if
       end do

       return
       end function maxt


C ===============================================================

