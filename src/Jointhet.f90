!==========================================================
!
!      Joint Latent class mixed model for continuous
!          Gaussian outcome
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



!----------------------------------------------------------
!
!- Module COMMUN avec les donnees dynamiques
!
!----------------------------------------------------------
    module commun_joint
      integer,save::ns
      double precision,dimension(:),allocatable,save::Y,time_cvpl
      double precision,dimension(:,:),allocatable,save::X
      integer,dimension(:),allocatable,save::nmes
      double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
      integer,dimension(:),allocatable,save::Devt
      integer,dimension(:),allocatable,save::ind_survint
    end module commun_joint
    
    
    
    module commun_modele_joint
      integer,save ::ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,logspecif &
           ,nprob,nrisq,nobs,nvarprob,maxmes,nmes_curr,nmes_curr_s
      integer,dimension(:),allocatable,save ::idxevt,idvdep
      double precision, dimension(:),allocatable,save::zi
      integer,dimension(:),allocatable,save::idea,idg,idprob
      integer,dimension(:),allocatable,save::prior
      integer ,save::risqcom,nvdepsurv,nvarxevt,typrisq,&
           nprisq,nz,idtrunc
      double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
           Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0 &
           ,Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
           Timt2,Timt3
      double precision,dimension(:),allocatable,save::Tmm_est,Tmm1_est,&
           Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est
      double precision,dimension(:),allocatable,save::Tmm_valt,Tmm1_valt, &
           Tmm2_valt,Tmm3_valt,Tim_valt,Tim1_valt,Tim2_valt,Tim3_valt
      integer,dimension(:),allocatable::indT
    end module commun_modele_joint
    
    





!===========================================================
!
      !SUBROUTINES
!
!===========================================================


!-----------------------------------------------------------
!                    SUBROUTINE PRINCIPALE: JOINTHET
!------------------------------------------------------------



      subroutine jointhet(Y0,X0,Prior0,idprob0,idea0,idg0,idxevt0 &
      ,ns0,ng0,nv0,nobs0   &
      ,nea0 &
        ,nmes0,idiag0,nwg0,npm0,b,vopt,vrais,ni     &
      ,istop,gconv,ppi0,ppitest0,resid_m0,resid_ss0,pred_m_g0     &
      ,pred_ss_g0   &
      ,pred_RE &
      ,convB,convL,convG,maxiter0,typrisq0        &
      ,idtrunc0,risqcom0,nz0,zi0,Tentr0,Tevt0          &
      ,Tsurvint0,devt0,ind_survint0,statsc0,risq_est,risqcum_est,nsim &
      ,time,logspecif0)




      use optim
      use commun_joint
      use commun_modele_joint
      use parameters

      IMPLICIT NONE


    !Declaration des variables en entree
      integer,intent(in):: nv0, maxiter0,nsim,nea0,logspecif0
      integer, intent(in) :: ns0, ng0, nobs0, idiag0, nwg0, npm0
      integer, dimension(nv0), intent(in) :: idea0,idg0,idprob0
      integer,dimension(nv0),intent(in) ::idxevt0
      integer, dimension(ns0), intent(in) :: nmes0,prior0
      double precision,dimension(ns0):: Tentr0,Tevt0,Tsurvint0
      double precision,dimension(nobs0),intent(in):: Y0
      double precision,dimension(nobs0*nv0), intent(in) :: X0
      double precision, intent(in) :: convB, convL, convG
      integer::risqcom0,typrisq0,nz0,idtrunc0
      integer,dimension(ns0) ::Devt0,ind_survint0
      double precision, dimension(nz0),intent(in)::zi0

    !Declaration des variable en entree et sortie
      double precision, dimension(npm0), intent(inout) :: b
    !Declaration des variables en sortie
      double precision, dimension(npm0*(npm0+1)/2),intent(out) :: vopt
      double precision,intent(out)::statsc0
      double precision :: vrais
      double precision, dimension(3), intent(out) :: gconv
      double precision, dimension(ns0*ng0), intent(out) :: ppi0
      double precision, dimension(ns0*ng0), intent(out) :: ppitest0
      double precision, dimension(nobs0), intent(out) :: resid_m0
      double precision, dimension(nobs0), intent(out) :: resid_ss0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_m_g0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_ss_g0
      double precision, dimension(nsim*ng0), intent(out) :: risq_est,risqcum_est
      integer, intent(inout) :: ni, istop
      double precision, dimension(nsim), intent(out)::time
      double precision, dimension(ns0*nea0), intent(out)::pred_RE



    !Variables locales
      integer :: jtemp,i,g,j,npm,ij,ier,k,ktemp,ig,nmestot,it
      double precision :: eps, ca, cb, dd
      double precision, dimension(npm0) :: mvc
      double precision, dimension(npm0*(npm0+3)/2) :: V
      double precision, dimension(nobs0) :: resid_m, resid_ss
      double precision, dimension(nobs0,ng0):: pred_m_g, pred_ss_g
      double precision, dimension(ns0,ng0)::ppi,ppitest
      double precision,dimension(:),allocatable::brisq_est
      integer::nef
      double precision,external::funcpaj


!      write(*,*)'entree',npm0,logspecif0


! *********************************
! initialisation des sorties
      resid_m0 = 0.d0
      pred_m_g0 = 0.d0
      resid_ss0 = 0.d0
      pred_ss_g0 = 0.d0
      ppi0 =  0.d0
      ppitest0 =  0.d0
      statsc0 =  9999.d0
      risq_est =  0.d0
      risqcum_est =  0.d0
      pred_RE=0.d0
      vopt=0.d0
      gconv=0.d0
      time(1)=zi0(1)
      do j=2,nsim
        time(j)=time(j-1)+(zi0(nz0)-zi0(1))/dble(nsim-1)
      end do
      time(nsim)=zi0(nz0)
! *********************************



      allocate(idprob(nv0)  &
     ,idea(nv0),idg(nv0),nmes(ns0),Tsurv0(ns0),Tsurv(ns0)    &
     ,Tsurvint(ns0),ind_survint(ns0),idxevt(nv0)             &
     ,idvdep(nv0),devt(ns0),prior(ns0))



      typrisq=typrisq0
      risqcom=risqcom0
      idtrunc=idtrunc0
      Tsurv0=Tentr0
      Tsurv=Tevt0
      Tsurvint=Tsurvint0
      devt=devt0
      ind_survint=ind_survint0
      logspecif=logspecif0

      maxiter=maxiter0
      maxmes=0

      nvdepsurv=0
      if(sum(ind_survint).gt.0) then
         nvdepsurv=1
      end if

      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do


      epsa=convB
      epsb=convL
      epsd=convG

!----------------- allocation generale ---------------------------


      allocate(Y(nobs0),X(nobs0,nv0))


!      nHT=0
!      do i=1,ns0
!         if(tsurv(i).gt.tsurvint(i).and.Devt(i).eq.1) then
!            nHT=nHT+1
!         end if
!
!         if(tsurv(i).gt.tsurvint(i).and.ind_survint(i).ne.1) then
!            write(*,*)'probleme HT',i,ind_survint(i),tsurvint(i),tsurv(i)
!         end if

!         if(tsurv(i).le.tsurvint(i).and.ind_survint(i).ne.0) then
!            write(*,*)'probleme HT',i,ind_survint(i),tsurvint(i),tsurv(i)
!         end if
!      end do

!      write(*,*)'sum',sum(ind_survint),nvdepsurv,nHT

!      nHT=0
!      do i=1,ns0
!         if(tsurv(i).gt.tsurvint(i).and.Devt(i).eq.0) then
!            nHT=nHT+1
!         end if
!      end do
!      write(*,*)'sum2',nHT

!------------------------
      
      nz=nz0
      select case (typrisq)
      case (1)
         allocate(zi(nz0))
         zi(1:nz0)=zi0(1:nz0)
      case (2)
         allocate(zi(nz0))
         zi(1:nz0)=zi0(1:nz0)
      case (3)
         allocate(zi(-2:nz0+3))
         do i=1,nz
            zi(i)=zi0(i)
         end do
      end select
      
      prior=0
      ppi=0.d0
      ppitest=0.d0
      resid_m0=0.d0
      resid_ss0=0.d0
      pred_m_g0=0.d0
      pred_ss_g0=0.d0
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
!------------------------
      if (typrisq.eq.1) then
         nprisq=nz-1
      end if
      if (typrisq.eq.2) then
         nprisq=2
      end if
      if (typrisq.eq.3) then
         nprisq=nz+2
      end if
!--------------------
      nmes=nmes0
      Y=0.d0
      X=0.d0
      idprob=0
      idea=0
      idg=0
      idxevt=0
      nmestot=0
      ktemp=0

      do k=1,nv
         idprob(k)=idprob0(k)
         idea(k)=idea0(k)
         idg(k)=idg0(k)
         idxevt(k)=idxevt0(k)
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
      
      nrisq=0
      
      if (risqcom.eq.1) then
         nrisq=nprisq
      end if
      if (risqcom.eq.2) then
         nrisq=nprisq+ng-1
      end if
      if (risqcom.eq.0) then
         nrisq=nprisq*ng
      end if
      ! parms fixes pour les vexp dans la modelisation de l'evenement
      
      nvarxevt=nvdepsurv
      Do k=1,nv
         If (idxevt(k).eq.1) then
            nvarxevt=nvarxevt+1
         end if
         If (idxevt(k).eq.2) then
            nvarxevt=nvarxevt+ng
         end if
      end do
      

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
         go to 1588
      end if


!  nb effets fixes = nb effets fixes sans melange
!                  + ng fois le nb de var dans melange

      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

      nef=nprob+nvarxevt+nrisq+ncssg+ncg*ng !
      npm=nef+nvc+nwg+1



!    write(*,*)'npm',npm,npm0
      if (npm.ne.npm0)then
         istop=4
         goto 1589
      end if
      
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
! Appel de la fonction spline
      if (typrisq.eq.3) then

      allocate(Tmm(ns0),Tmm1(ns0),Tmm2(ns0),Tmm3(ns0),Tim(ns0)         &
     ,Tim1(ns0),Tim2(ns0),Tim3(ns0),Tmm0(ns0),Tmm01(ns0),Tmm02(ns0)  &
     ,Tmm03(ns0),Tim0(ns0),Tim01(ns0),Tim02(ns0),Tim03(ns0),          &
      Tmmt(ns0),Tmmt1(ns0),Tmmt2(ns0),Tmmt3(ns0),Timt(ns0),Timt1(ns0) &
     ,Timt2(ns0),Timt3(ns0))

          call splines()
      end if

! lancement de l'optimisation

      ca=0.d0
      cb=0.d0
      dd=0.d0
!      write(*,*)'B avant marq98',npm,(b(i),i=1,npm)
      call marq98(b,npm,ni,v,vrais,ier,istop,ca,cb,dd,funcpaj)




!        write(*,*)'B apres marq98',(b(i),i=1,npm)
!           write(*,*)'-------fin de optimisation avec marq98-------'
!               write(*,*)'ier',ier,'istop',istop
!           write(*,*)'ca',ca,'cb',cb,'dd',dd
!
!               write(*,*)
!               write(*,*)'    FIN OPTIMISATION  ..... '
!               write(*,*)'istop',istop,'vrais',vrais





      gconv=0.d0
      gconv(1)=ca
      gconv(2)=cb
      gconv(3)=dd

      vopt=0.d0
      vopt(1:npm0*(npm0+1)/2)=v(1:npm0*(npm0+1)/2)

      do k=1,nwg+1
         b(nef+nvc+k)=abs(b(nef+nvc+k))
      end do

      ! resultats : SE parms transformes


      if (istop==4 .or. istop==12)then
         goto 1589
      end if

      if (istop.eq.1.or.istop.eq.3) then
         !------ estimation de la variance des probabilites des composantes

         if (ng.gt.1) then
            call postprobj(B,npm,ppi,ppitest)
         end if
         
         call residualsj(b,npm,ppi,resid_m,pred_m_g,resid_ss &
              ,pred_ss_g,pred_RE)
         
         
         ig=0
         ij=0
         do i=1,ns
            do g=1,ng0
               ig=ig+1
               ppi0(ig)=PPI(i,g)
               ppitest0(ig)=ppitest(i,g)
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

         call scoretest(b,npm,statsc0)

         if (typrisq.eq.3) then

            allocate(Tmm_est(nsim),Tmm1_est(nsim),Tmm2_est(nsim),Tmm3_est(nsim)&
                 ,Tim_est(nsim),Tim1_est(nsim),Tim2_est(nsim),Tim3_est(nsim))
            call splines_estime(time,nsim)
         end if

         allocate(brisq_est(nprisq))

         do g=1,ng
            brisq_est=0.d0
            if (logspecif.eq.1.or.typrisq.eq.2) then
                  if (risqcom.eq.0) then
               do k=1,nprisq
                    brisq_est(k)=exp(b(nprob+nprisq*(g-1)+k))
                  end do
                  else if (risqcom.eq.1) then
               do k=1,nprisq
                    brisq_est(k)=exp(b(nprob+k))
                  end do
                  else if (risqcom.eq.2) then
               do k=1,nprisq
                    brisq_est(k)=exp(b(nprob+k))
                  end do
                  end if
            else
                  if (risqcom.eq.0) then
               do k=1,nprisq
                     brisq_est(k)=b(nprob+nprisq*(g-1)+k) &
                          *b(nprob+nprisq*(g-1)+k)
                  end do
                  elseif (risqcom.eq.1) then
               do k=1,nprisq
                     brisq_est(k)=b(nprob+k)  &
                          *b(nprob+k)
                  end do
                  elseif (risqcom.eq.2) then
               do k=1,nprisq
                     brisq_est(k)=b(nprob+k)  &
                          *b(nprob+k)
                  end do
                end if
            end if
            call fct_risq_estime(brisq_est,time,nsim,g,risq_est,risqcum_est)

            if (risqcom.eq.2.and.ng.gt.1.and.g.lt.ng) then
               do i=1,nsim
                  risq_est(nsim*(g-1)+i)=risq_est(nsim*(g-1)+i)*exp(b(nprob+nprisq+g))
                  risqcum_est(nsim*(g-1)+i)=risqcum_est(nsim*(g-1)+i)*exp(b(nprob+nprisq+g))
               end do
            end if
         end do

         deallocate(brisq_est)

         if (typrisq.eq.3) then
            deallocate(Tmm_est,Tmm1_est,Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est)
         end if

      end if


 1589 continue


      if (typrisq.eq.3) then
      deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
      Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
      Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
      endif

 1588 continue

!--------------- liberation de la memoire ------------------------
      deallocate(Y,X,idprob,idea,idg,nmes,Tsurv0,Tsurv,Tsurvint &
     , ind_survint,zi,idxevt,idvdep,devt,prior)


! 100  continue


      return


      end subroutine jointhet

!-----------------------------------------------------------
!                        FUNCPA
!------------------------------------------------------------

      double precision function funcpaj(b,m,id,thi,jd,thj)
        
        
        use commun_joint,only:ns,nmes
        use commun_modele_joint,only:nmes_curr
        implicit none
        
        integer::m,i,id,jd
        double precision::rl,thi,thj,funcpij,rli
        double precision,dimension(m)::b
        
        !    write(*,*)'entree funcpa'
        
        rl=0.d0
        rli=0.d0
        nmes_curr=0
        do i=1,ns
           rli =  funcpij(b,m,id,thi,jd,thj,i)
           if (rli.eq.-1.d9) then
              rl = -1.d9
              goto 468
           end if
           rl = rl + rli
           nmes_curr = nmes_curr + nmes(i)
           
        end do
468     continue
        funcpaj =rl
        !    write(*,*)'sortie funcpa',funcpaj
        
        return
        
      end function funcpaj
      






!-----------------------------------------------------------
!                        FUNCPI
!------------------------------------------------------------

      
      double precision function funcpij(b,npm,id,thi,jd,thj,i)
        use commun_joint
        use commun_modele_joint
        use optim
        
        IMPLICIT NONE
        
        integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm
        integer ::ier,nmoins
        double precision,dimension(maxmes,nv) ::Z,P,X00,X2
        double precision,dimension(nvarprob) ::Xprob,bprob
        double precision,dimension(nv,nv) ::Ut,Ut1
        !      double precision,dimension(maxmes,maxmes) ::VC,Se
        double precision,dimension(:,:),allocatable ::VC,Se
        double precision,dimension(npm) :: b,b1
        !      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
        double precision,dimension(:),allocatable ::Vi
        double precision,dimension(nv) :: b0,b2
        double precision :: vrais,eps,det
        double precision ::thi,thj,temp
        double precision ::Y4,expo
        double precision,dimension(maxmes) :: mu,Y1,Y2,Y3
        double precision,dimension(ng) :: pi
        double precision,dimension(nrisq)::brisq
        double precision::retard,entretard,vrais_survie
        double precision,dimension(nvarxevt)::bevt,Xevt
        double precision::bevtint
        double precision,dimension(ng)::risq,surv,surv0,  &
             survint
        logical::is_nan
        integer::nef
        
        !    write(*,*)'entree funcpi',i
        
        
        allocate(VC(maxmes,maxmes),Se(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
        VC=0.d0
        Vi=0.d0
        Se=0.d0
        funcpij=0.d0
        
        
        b1=0.d0
        eps=1.D-20
        do k=1,npm
           b1(k)=b(k)
        end do
        if (id.ne.0) b1(id)=b1(id)+thi
        if (jd.ne.0) b1(jd)=b1(jd)+thj
        
        !----------- rappel des parametres utilises ---------
        
        nef=nprob+nrisq+nvarxevt+ncssg+ncg*ng
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
        
        ! calcul de brisq e chaque composante et risq, surv et surv0 pour tous les i et tous les g
        
        risq=0.d0
        surv=0.d0
        surv0=0.d0
        survint=0.d0
        
        do g=1,ng
           
           brisq=0.d0
           if (logspecif.eq.1.or.typrisq.eq.2) then
              if (risqcom.eq.0) then
                 do k=1,nprisq
                    brisq(k)=exp(b1(nprob+nprisq*(g-1)+k))
                 end do
              elseif(risqcom.eq.1) then
                 do k=1,nprisq
                    brisq(k)=exp(b1(nprob+k))
                 end do
              elseif (risqcom.eq.2) then
                 do k=1,nprisq
                    brisq(k)=exp(b1(nprob+k))
                 end do
              end if
              
           else
              if (risqcom.eq.0) then
                 do k=1,nprisq
                    brisq(k)=b1(nprob+nprisq*(g-1)+k) &
                         *b1(nprob+nprisq*(g-1)+k)
                 end do
              elseif(risqcom.eq.1) then
                 do k=1,nprisq
                    brisq(k)=b1(nprob+k)*b1(nprob+k)
                 end do
              elseif (risqcom.eq.2) then
                 do k=1,nprisq
                    brisq(k)=b1(nprob+k)*b1(nprob+k)
                 end do
              end if
           end if
           
           
           call fct_risq_i(i,brisq,g,risq,surv,surv0,survint)
           
           
           if (risqcom.eq.2.and.ng.gt.1.and.g.lt.ng) then
              risq(g)=risq(g)*exp(b1(nprob+nprisq+g))
              surv(g)=surv(g)*exp(b1(nprob+nprisq+g))
              survint(g)=survint(g)*exp(b1(nprob+nprisq+g))
              surv0(g)=surv0(g)*exp(b1(nprob+nprisq+g))
           end if
        end do
        
        !      if(tsurvint(i).ne.tsurv(i).and.ind_survint(i).eq.0) then
        !         write(*,*)'problem',tsurvint(i),tsurv(i)
        !      end if
        !      do g=1,ng
        !         if(surv(g).ne.survint(g).and.ind_survint(i).eq.0) then
        !            write(*,*)'problem surv',surv(g),survint(g)
        !         end if
        !         if(surv(g).eq.survint(g).and.ind_survint(i).eq.1) then
        !            write(*,*)'problem surv ind',surv(g),survint(g)
        !         end if
        !      end do
        
        
        
        
        
        ! ----------- boucle sur les individus -------------
        entretard=0.d0
        vrais_survie=0.d0
        vrais=0.d0
        
        ! -------- creation de Vi = ZiGZi'+se*seIni ----------
        ! creation de Zi
        
        Z=0.d0
        l=0
        do k=1,nv
           if (idea(k).eq.1) then
              l=l+1
              do j=1,nmes(i)
                 Z(j,l)=dble(X(nmes_curr+j,k))
              end do
           end if
           
        end do
        ! creation de s2*I et Y1
        
        Se=0.d0
        Y1=0.d0
        do j=1,nmes(i)
           !      kk=kk+1
           Se(j,j)=b1(npm)*b1(npm)
           Y1(j)=dble(Y(nmes_curr + j))
        end do
        
        ! creer Xevt:
        Xevt=0.d0
        
        l=0
        do k=1,nv
           if (idxevt(k).eq.1.or.idxevt(k).eq.2) then
              l=l+1
              Xevt(l)=X(nmes_curr+1,k)
           end if
        end do
        
        ! creation de P=Zi*Ut et V=P*P' que si non spec aux classes
        
        IF (nwg.eq.0.OR.NG.EQ.1) then
           
           
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
           CALL DSINV(Vi,nmes(i),eps,ier,det)
           if (ier.eq.-1) then
              !               write(*,*)'probleme dsinv'
              funcpij=-1.d9
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
        END IF
        
        !     debut du calcul de la vraisemblance
        vrais=vrais-nmes(i)*dlog(dble(2*3.14159265))
        
        
        
! Cas 1 : ng=1
        IF (ng.eq.1) then
           
           bevt=0.d0
           
           bevtint=0.d0
           l=1
           do k=1,nv!warning nvar
              if (idxevt(k).eq.1) then
                 bevt(l)=b1(nprob+nrisq+l)
                 l=l+1
              end if
           end do
           
           if (l-1.ne.nvarxevt-nvdepsurv) then
              !              write(*,*)'probleme nvarxevt'
              funcpij=-1.d9
              goto 654
           end if
           
           if (ind_survint(i).eq.1) then
              bevtint=b1(nprob+nrisq+nvarxevt)
           end if
           
           
           !--------------------------------------------------------------
           
           vrais=vrais-det
           b0=0.d0
           l=0
           X00=0.d0
           do k=1,nv
              if (idg(k).ne.0) then
                 l=l+1
                 do j=1,nmes(i)
                    X00(j,l)=dble(X(nmes_curr+j,k))
                 end do
                 b0(l)=b1(nprob+nrisq+nvarxevt+l)
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
           
           !------------------------------------------------------------------


           if (Devt(i).eq.1) then
              if (ind_survint(i).eq.1) then
                 risq(1)=risq(1)*exp(bevtint)
              end if
              vrais=vrais+2*(log(risq(1))+DOT_PRODUCT(Xevt,bevt) &
                   -exp(DOT_PRODUCT(Xevt,bevt))*(survint(1)+     &
                   exp(bevtint)*(surv(1)-survint(1))))
              
              
              vrais_survie=vrais_survie+2*(log(risq(1))        &
                   +DOT_PRODUCT(Xevt,bevt)                        &
                   -exp(DOT_PRODUCT(Xevt,bevt))*(survint(1)+  &
                   exp(bevtint)*(surv(1)-survint(1))))
              !        else
              !           vrais=vrais+2*(log(risq(1))+DOT_PRODUCT(Xevt,bevt) &
              !                -exp(DOT_PRODUCT(Xevt,bevt))*surv(1))
              !
              !           vrais_survie=vrais_survie+2*(log(risq(1))        &
              !                +DOT_PRODUCT(Xevt,bevt)                        &
              !                -exp(DOT_PRODUCT(Xevt,bevt))*(surv(1)))
              !        end if
              if (risq(1).le.0.or.is_nan(risq(1))) then
                 !               write(*,*)'probleme risq',risq(i,1)
                 funcpij=-1.d9
                 goto 654
              end if
              
              !           if(ind_survint(i).eq.1.and.thi.eq.0.and.thj.eq.0) then
              !              write(*,*)'devt1',i,surv(1),survint(1),tsurv(i),Tsurvint(i)
              !           end if
              
           end if
           if (Devt(i).eq.0) then
              !           if (ind_survint(i).eq.1) then
              vrais=vrais-2*exp(DOT_PRODUCT(Xevt,bevt)) &
                   *(survint(1)+exp(bevtint)       &
                   *(surv(1)-survint(1)))
              
              vrais_survie=vrais_survie                 &
                   -2*exp(DOT_PRODUCT(Xevt,bevt))      &
                   *(survint(1)+exp(bevtint)      &
                   *(surv(1)-survint(1)))
              !           else
!              vrais=vrais-2*exp(DOT_PRODUCT(Xevt,bevt)) &
!                   *(surv(1))
!
!              vrais_survie=vrais_survie                 &
!                   -2*exp(DOT_PRODUCT(Xevt,bevt))      &
!                   *(surv(1))
 !          end if

!           if(ind_survint(i).eq.1.and.thi.eq.0.and.thj.eq.0) then
!              write(*,*)'devt0',i,surv(1),survint(1),tsurv(i),Tsurvint(i)
!           end if

           end if
           entretard=entretard-surv0(1)*exp(DOT_PRODUCT(Xevt,bevt))
           
           
           
           
        ELSE
           


! cas 2 :  ng>1  composantes
! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))

           if (prior(i).ne.0) then
              pi=0.d0
              pi(prior(i))=1.d0
           else
              Xprob=0.d0
              Xprob(1)=1
              l=0
              do k=1,nv
                 if (idprob(k).eq.1) then
                    l=l+1
                    Xprob(1+l)=X(nmes_curr+1,k)
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
           retard=0.d0
           ! creation des vecteurs de variables explicatives
           l=0
           m=0
           X00=0.d0
           X2=0.d0
           do k=1,nv
              if (idg(k).eq.2) then
                 l=l+1
                 do j=1,nmes(i)
                    X2(j,l)=dble(X(nmes_curr+j,k))
                 end do
              else if (idg(k).eq.1) then
                 m=m+1
                 do j=1,nmes(i)
                    X00(j,m)=dble(X(nmes_curr+j,k))
                 end do
              end if
           end do
           
           b2=0.d0
           b0=0.d0
           expo=0.d0
           
           
           Do g=1,ng
              !-------------------------------------------------------
              
              ! bevt
              bevt=0.d0
              bevtint=0.d0
              
              m=0
              l=1
              do k=1,nv
                 
                 if (idxevt(k).eq.1) then
                    
                    bevt(l)=b1(nprob+nrisq+m+1)
                    l=l+1
                    m=m+1
                 end if
                 if (idxevt(k).eq.2) then
                    bevt(l)=b1(nprob+nrisq+m+g)
                    l=l+1
                    m=m+ng
                 end if
                 
              end do
              if (nvdepsurv.ne.0) then
                 bevtint=b1(nprob+nrisq+nvarxevt)
              end if
              
!     write(*,*)'g=',g,'bevt=',(bevt(k),k=1,l)
!-------------------------------------------------------

               nmoins=0
               l2=0
               m2=0
               do k=1,nv
                  if (idg(k).eq.1) then
                     m2=m2+1
                     b0(m2)=b1(nprob+nrisq+nvarxevt+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nrisq+nvarxevt+nmoins+g)
                     nmoins=nmoins+ng
                  end if
               end do

! variance covariance si spec aux classes :
               Ut1=0.d0
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

                  CALL DSINV(Vi,nmes(i),eps,ier,det)
                  if (ier.eq.-1) then
!               write(*,*)'probleme dsinv'
                     funcpij=-1.d9
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
               Y4=DOT_PRODUCT(Y2,Y3)


!           if (i.lt.3.and.thi.eq.0.and.thj.eq.0) then
!              write(*,*)'i',i,g,devt(i),ind_survint(i)
!              write(*,*)'i',i,g,risq(i,g),survint(i,g),surv(i,g)
!              write(*,*)'i',i,g,tsurv(i),tsurv0(i),tsurvint(i)
!             write(*,*)'bevt',g,bevt
!           end if


                  if (ind_survint(i).eq.1) then
                     risq(g)=risq(g)*exp(bevtint)
                  end if

                  if (Devt(i).eq.1) then
                     expo=expo+pi(g)*risq(g)*                        &
                          exp(DOT_PRODUCT(Xevt,bevt)                    &
                          +(-det-Y4)/2.d0-exp(DOT_PRODUCT(Xevt,bevt))*  &
                          (survint(g)+exp(bevtint)*(surv(g)      &
                          -survint(g))))
                  end if
                  if (Devt(i).eq.0) then
                     expo=expo+pi(g)*exp((-det-Y4)/2.d0                &
                          -exp(DOT_PRODUCT(Xevt,bevt))*                 &
                          (survint(g)+exp(bevtint)*(surv(g)      &
                          -survint(g))))
                  end if


               retard=retard+pi(g)*exp(-surv0(g) &
                    *exp(DOT_PRODUCT(Xevt,bevt)))

            end do



            if (expo.le.0.or.is_nan(expo).or.expo.gt.1.d30) then
!               write(*,*)'i',i,Devt(i)
!               write(*,*)'Y',Y4,det
!               write(*,*)'Y',nmes(i),(Y1(j),j=1,nmes(i))
!               do l=1,nea
!                  write(*,*)'Z',nmes(i),(Z(j,l),j=1,nmes(i))
!               end do
!               write(*,*)'expo le 0',expo,(pi(g),g=1,ng)
!               write(*,*)'risq(i,g)',(risq(i,g),g=1,ng)
!               write(*,*)'surv(i,g)',(surv(i,g),g=1,ng)
!               write(*,*)'survint(i,g)',(survint(i,g),g=1,ng)
!               write(*,*)'expo lt 0 funcpa',expo
               funcpij=-1.d9
               goto 654
            end if

            entretard=entretard+log(retard)

            vrais=vrais+2*log(expo)

         END IF

      if (idtrunc.eq.0) then
         entretard=0.d0
      end if


      if (ng.eq.1) then
         vrais_survie=vrais_survie/2.d0-entretard
      end if

      if (vrais.lt.-1.d9) then
!         write(*,*)'vrais inf -10^9',vrais
         funcpij=-1.d9
         goto 654
       end if


! FIN BOUCLE SUJET
      funcpij=vrais/2.D0-entretard

 654  continue

      deallocate(VC,Se,Vi)


      if (is_nan(funcpij).or.abs(funcpij).gt.1.d30) then
!         write(*,*)'is_nan funcpa',vrais
          funcpij=-1.d9
       end if
      return

      end function funcpij





!------------------------------------------------------------
!                      POSTPROB
!------------------------------------------------------------

!-------------------------------------------------------------
!
!          Subroutine pour calculer les
!      probabilites a posteriori de suivre chacune
!      des composantes g pour chacun des sujets i
!-------------------------------------------------------------

      subroutine postprobj(b,npm,ppi,ppiy)

    use commun_joint
    use commun_modele_joint
      use optim

      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk,nef
      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision ::Y4,fevt
      double precision,dimension(ng) ::fi,pi,fi1
      double precision,dimension(ns,ng) ::ppiy,ppi
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3

      double precision,dimension(nvarxevt)::bevt,Xevt
      double precision,dimension(ns,ng)::risq2
      double precision::bevtint
      double precision,dimension(nprisq)::brisq
      double precision,dimension(ns,ng)::risq,surv,surv0,  &
      survint

      eps=1.D-20
      ppiy=0.D0
      ppi=0.D0
      b1=0.d0
      do k=1,npm
         b1(k)=b(k)
      end do

      allocate(VC(maxmes,maxmes),Se(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      Vi=0.d0
      Se=0.d0
      risq2=0.d0



!----------- rappel des parametres utilises ---------

      nef=nprob+nrisq+nvarxevt+ncssg+ncg*ng

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

!------------------- calcul de risq ------------------------
! calcul de brisq e chaque composante et risq, surv et surv0 pour tous les i et tous les g


      risq=0.d0
      surv=0.d0
      surv0=0.d0
      survint=0.d0

      do g=1,ng

         brisq=0.d0
         if (logspecif.eq.1.or.typrisq.eq.2) then
               if (risqcom.eq.0) then
                 do k=1,nprisq
                   brisq(k)=exp(b1(nprob+nprisq*(g-1)+k))
               end do
               elseif (risqcom.eq.1) then
                do k=1,nprisq
                  brisq(k)=exp(b1(nprob+k))
                end do
               elseif (risqcom.eq.2) then
             do k=1,nprisq
                brisq(k)=exp(b1(nprob+k))
              end do
            endif
        else
               if (risqcom.eq.0) then
              do k=1,nprisq
                  brisq(k)=b1(nprob+nprisq*(g-1)+k)   &
                      *b1(nprob+nprisq*(g-1)+k)
               end do
               elseif (risqcom.eq.1) then
                do k=1,nprisq
                  brisq(k)=b1(nprob+k) *b1(nprob+k)

                end do
               elseif (risqcom.eq.2) then
               do k=1,nprisq
                brisq(k)=b1(nprob+k)*b1(nprob+k)
               end do
              end if
        end if
!         write(*,*)'g=',g,'brisq=',(brisq(k),k=1,nrisq)

         call fct_risq(brisq,g,risq,surv,surv0,survint)

         if (risqcom.eq.2.and.ng.gt.1.and.g.lt.ng) then
            do i=1,ns
               risq(i,g)=risq(i,g)*exp(b1(nprob+nprisq+g))
               surv(i,g)=surv(i,g)*exp(b1(nprob+nprisq+g))
               survint(i,g)=survint(i,g)*exp(b1(nprob+nprisq+g))
               surv0(i,g)=surv0(i,g)*exp(b1(nprob+nprisq+g))
            end do
         end if

      end do

! ----------- boucle sur les individus -------------
      kk=0
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


         Xevt=0.d0
            l=0
            do k=1,nv
               if (idxevt(k).eq.1.or.idxevt(k).eq.2) then
                  l=l+1
                  Xevt(l)=X(it+1,k)
               end if

            end do
! creation de s2*I et Y1

         Se=0.d0
     Y1=0.d0
         do j=1,nmes(i)
        kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
        Y1(j)=dble(Y(kk))
         end do

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

         CALL DSINV(Vi,nmes(i),eps,ier,det)
         if (ier.eq.-1) then
            ppiy=-1.d0
        ppi=-1.d0
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

! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       if (prior(i).ne.0) then
          pi=0.d0
          pi(prior(i))=1.d0
       else
           Xprob=0.d0
           Xprob(1)=1
           l=0
           do k=1,nv
             if (idprob(k).eq.1) then
                l=l+1
                Xprob(1+l)=X(it+1,k)
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
 
       fi=0.d0
       fi1=0.d0
       b0=0.d0
       b2=0.d0

       do g=1,ng
! bevt
            bevt=0.d0
            bevtint=0.d0
               m=0
               l=1
               do k=1,nv
                  if (idxevt(k).eq.1) then
                     bevt(l)=b1(nprob+nrisq+m+1)
                     l=l+1
                     m=m+1
                  end if
                  if (idxevt(k).eq.2) then
                     bevt(l)=b1(nprob+nrisq+m+g)
                     l=l+1
                     m=m+ng
                  end if

               end do
               if (nvdepsurv.ne.0) then
                  bevtint=b1(nprob+nrisq+nvarxevt)
               end if

          nmoins=0
          l2=0
          m2=0
          do k=1,nv
             if (idg(k).eq.1) then
                m2=m2+1
                b0(m2)=b1(nprob+nrisq+nvarxevt+nmoins+1)
                nmoins=nmoins+1
             else if (idg(k).eq.2) then
                l2=l2+1
                b2(l2)=b1(nprob+nrisq+nvarxevt+nmoins+g)
                nmoins=nmoins+ng
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

             CALL DSINV(Vi,nmes(i),eps,ier,det)
             if (ier.eq.-1) then
                ppi=-1.d0
        ppiy=-1.d0
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

          risq2(i,g)=risq(i,g)
          if (ind_survint(i).eq.1) then
             risq2(i,g)=risq(i,g)*exp(bevtint)
          end if

          fi(g)=fi(g)- nmes(i)*log(dble(2*3.14159265))
          fi(g)=fi(g) -det
          fi(g)=fi(g) - Y4
          fi(g)=fi(g)/(2.d0)
          fi1(g)=dexp(fi(g))


          fevt=log(risq2(i,g))+DOT_PRODUCT(bevt,Xevt)

              fi(g)=fi(g)+Devt(i)*fevt-exp(DOT_PRODUCT(Xevt,bevt))* &
              (survint(i,g)+exp(bevtint)*(surv(i,g)-survint(i,g)))

      fi(g)=exp(fi(g))


       end do
       do g=1,ng
          ppiy(i,g)=pi(g)*fi1(g)/DOT_PRODUCT(pi,fi1)
      ppi(i,g)=pi(g)*fi(g)/DOT_PRODUCT(pi,fi)
       end do

       it=it+nmes(i)

      end do


 147  continue


      deallocate(VC,Vi,Se)

      return

      end subroutine postprobj






!------------------------------------------------------------
!                      RESIDUALS
!------------------------------------------------------------


      subroutine residualsj(b1,npm,ppi,resid_m,pred_m_g,resid_ss, &
      pred_ss_g,pred_RE)
        use commun_joint
        use commun_modele_joint
        use optim
        implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,npm
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk
      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob,err2
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Se,VC1
      double precision,dimension(npm) ::b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
!      double precision,dimension(nv) ::b3
      double precision :: eps,det
      double precision ::temp
      double precision,dimension(nv,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m &
      ,pred_m,resid_ss,pred_ss
      double precision,dimension(ns*nea)::pred_RE
      double precision,dimension(nobs,ng)::pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
    integer::nef



      allocate(VC(maxmes,maxmes),VC1(maxmes,maxmes),Se(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      VC1=0.d0
      Vi=0.d0
      Se=0.d0



      eps=1.D-20

!----------- rappel des parametres utilises ---------

      nef=nprob+nrisq+nvarxevt+ncssg+ncg*ng



!      write(*,*)'nvc',nvc,'nea',nea,'nwg',nwg,'nef',nef
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

!         write(*,*)'Ut',Ut

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m=0.d0
      pred_m_g=0.d0
      pred_ss=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0
      pred_RE=0.d0

      nmes_cur=0
      kk=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+se*seIni ----------

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

!     creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0
         do j=1,nmes(i)
            kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
            Y1(j)=dble(Y(kk))
         end do

!     creation de P=Zi*Ut et V=P*P' que si non spec aux classes
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
                  b0(l)=b1(nprob+nrisq+nvarxevt+l)
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



!     cas 2 :  ng>1  composantes
         else

!     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            if (prior(i).ne.0) then
               pi=0.d0
               pi(prior(i))=1.d0
            else
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
!               b3=0.d0
               nmoins2=0
               do k=1,nv
                  if (idg(k).eq.1) then
                     m2=m2+1
                     b0(m2)=b1(nprob+nrisq+nvarxevt+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nrisq+nvarxevt+nmoins+g)
                     nmoins=nmoins+ng
                  end if
                  IF (IDEA(k).EQ.1.and.idg(k).eq.1) THEN
                     n2=n2+1
!                     b3(n2)=b1(nprob+nrisq+nvarxevt+nmoins2+1)
                     nmoins2=nmoins2+1
                  else if(IDEA(k).EQ.1.and.idg(k).eq.2) THEN
                     n2=n2+1
!                     b3(n2)=b1(nprob+nrisq+nvarxevt+nmoins2+g)
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


!     Vi en vecteur

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
                  pred_m_g(nmes_cur+j,g)=mu(j)
                  pred_ss_g(nmes_cur+j,g)=pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j)   &
                      +ppi(i,g)*(Y1(j)-pred1(j))
                  pred_ss(nmes_cur+j)=pred_ss(nmes_cur+j)+ppi(i,g) &
                      *pred1(j)
                  pred_m(nmes_cur+j)=pred_m(nmes_cur+j)+pi(g)*mu(j)
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))

               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*err2(k)
               end do

            end do

         end if


         nmes_cur=nmes_cur+nmes(i)
      end do

 654  continue

      deallocate(VC,VC1,Vi,Se)
      return

      end subroutine residualsj




!---------------------------------------------------------------
!
      !FCT_RISQ
!
!---------------------------------------------------------------






      subroutine fct_risq(brisq,g,risq,surv,surv0,survint)
        use commun_joint
        use commun_modele_joint

      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


      double precision, dimension(nprisq)::brisq
      integer ::i,g,ll,kk,ii,l,j
      double precision,dimension(ns,ng):: risq,surv,surv0,survint
      double precision::som



      do i=1,ns

         if (typrisq.eq.2) then

            surv(i,g)=brisq(1)*(tsurv(i))**brisq(2)

            risq(i,g)=brisq(1)*brisq(2) &
                *(Tsurv(i))**(brisq(2)-1)
            if (idtrunc.eq.1) then
               surv0(i,g)=brisq(1)*(tsurv0(i))**brisq(2)
            end if
            if (ind_survint(i).eq.1) then
               survint(i,g)=brisq(1)*(tsurvint(i)) &
      **brisq(2)
            else
               survint(i,g)=surv(i,g)
            end if

         end if

         if (typrisq.eq.1) then
            do j=1,nz-1
               som=0.d0
               do l=1,j-1
                  som=som+brisq(l)*(zi(l+1)-zi(l))
               end do
               if (idtrunc.eq.1) then
                  if (Tsurv0(i).ge.zi(j).and.Tsurv0(i).le.zi(j+1)) then
                     surv0(i,g)=som+brisq(j)*(Tsurv0(i)-zi(j))
                  end if
               end if
               if (Tsurv(i).ge.zi(j).and.Tsurv(i).le.zi(j+1)) then
                  surv(i,g)=som+brisq(j)*(Tsurv(i)-zi(j))
                  risq(i,g)=brisq(j)
               end if
               if (ind_survint(i).eq.1) then
                  if (Tsurvint(i).ge.zi(j).and.Tsurvint(i).le.zi(j+1)) &
                       then
                     survint(i,g)=som+brisq(j)*(Tsurvint(i)-zi(j))
                  end if
               end if
            end do
            if (ind_survint(i).eq.0) then
               survint(i,g)=surv(i,g)
            end if
         end if

         if (typrisq.eq.3) then
            !------------ survie et risq pour Tsurv ----------------
            ll=0
            if (Tsurv(i).eq.zi(nz)) then
               ll=nz-1
            end if
            som=0.d0
            do kk=2,nz
               if ((Tsurv(i).ge.zi(kk-1)).and.(Tsurv(i).lt.zi(kk))) &
                    then
                  ll=kk-1
               end if
            end do
            if (ll.gt.1) then
               do ii=1,ll-1
                  som=som+brisq(ii)
               end do
            end if

            surv(i,g)=som+brisq(ll)*Tim3(i)+brisq(ll+1)*Tim2(i) &
                 +brisq(ll+2)*Tim1(i)+brisq(ll+3)*Tim(i)
            risq(i,g)=brisq(ll)*Tmm3(i)+brisq(ll+1)*Tmm2(i)     &
                 +brisq(ll+2)*Tmm1(i)+brisq(ll+3)*Tmm(i)

            !------------ survie et risq pour Tsurv0 ----------------

            if (idtrunc.eq.1) then
               ll=0
               if (Tsurv0(i).eq.zi(nz)) then
                  ll=nz-1
               end if
               som=0.d0
               do kk=2,nz
                  if ((Tsurv0(i).ge.zi(kk-1)).and.(Tsurv0(i).lt.zi(kk))) &
                       then
                     ll=kk-1
                  end if
               end do
!               if (ll.lt.1.or.ll.gt.nz-1) then
!                  write(*,*) 'probleme dans fct_risq splines'
!                  write(*,*) 'll=',ll,'T=',Tsurv0(i)
!                  stop
!               end if
               if (ll.gt.1) then
                  do ii=1,ll-1
                     som=som+brisq(ii)
                  end do
               end if

               surv0(i,g)=som+brisq(ll)*Tim03(i)+brisq(ll+1)*Tim02(i) &
                    +brisq(ll+2)*Tim01(i)+brisq(ll+3)*Tim0(i)

            end if

            !------------ survie et risq pour Tsurvint ----------------


            if (ind_survint(i).eq.1) then
               ll=0
               if (Tsurvint(i).eq.zi(nz)) then
                  ll=nz-1
               end if
               som=0.d0
               do kk=2,nz
                  if((Tsurvint(i).ge.zi(kk-1)).and.(Tsurvint(i).lt.zi(kk))) &
                       then
                     ll=kk-1
                  end if
               end do
!               if (ll.lt.1.or.ll.gt.nz-1) then
!                  write(*,*) 'probleme dans fct_risq splines'
!                  write(*,*) 'll=',ll,'T=',Tsurvint(i)
!                  stop
!               end if
               if (ll.gt.1) then
                  do ii=1,ll-1
                     som=som+brisq(ii)
                  end do
               end if

               survint(i,g)=som+brisq(ll)*Timt3(i)+brisq(ll+1)*Timt2(i) &
                    +brisq(ll+2)*Timt1(i)+brisq(ll+3)*Timt(i)

            else
               survint(i,g)=surv(i,g)
            end if

         end if
      end do


      end subroutine fct_risq

!---------------------------------------------------------------
!
      !FCT_RISQI
!
!---------------------------------------------------------------


      subroutine fct_risq_i(i,brisq,g,risq,surv,surv0,survint)
    use commun_joint
    use commun_modele_joint

      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


      double precision, dimension(nprisq)::brisq
      integer ::i,g,ll,kk,ii,l,j
      double precision,dimension(ng)::risq,surv,surv0,survint
      double precision::som


         if (typrisq.eq.2) then

            surv(g)=brisq(1)*(tsurv(i))**brisq(2)

            risq(g)=brisq(1)*brisq(2)*(tsurv(i))**(brisq(2)-1)
            if (idtrunc.eq.1) then
               surv0(g)=brisq(1)*(tsurv0(i))**brisq(2)
            end if
            if (ind_survint(i).eq.1) then
               survint(g)=brisq(1)*(tsurvint(i))**brisq(2)
            else
               survint(g)=surv(g)
            end if

         end if

         if (typrisq.eq.1) then
            do j=1,nz-1
               som=0.d0
               do l=1,j-1
                  som=som+brisq(l)*(zi(l+1)-zi(l))
               end do
               if (idtrunc.eq.1) then
                  if (Tsurv0(i).ge.zi(j).and.Tsurv0(i).le.zi(j+1)) then
                     surv0(g)=som+brisq(j)*(Tsurv0(i)-zi(j))
                  end if
               end if
               if (Tsurv(i).ge.zi(j).and.Tsurv(i).le.zi(j+1)) then
                  surv(g)=som+brisq(j)*(Tsurv(i)-zi(j))
                  risq(g)=brisq(j)
               end if
               if (ind_survint(i).eq.1) then
                  if (Tsurvint(i).ge.zi(j).and.Tsurvint(i).le.zi(j+1)) &
                       then
                     survint(g)=som+brisq(j)*(Tsurvint(i)-zi(j))
                  end if
               end if
            end do
            if (ind_survint(i).eq.0) then
               survint(g)=surv(g)
            end if
         end if

         if (typrisq.eq.3) then
            !------------ survie et risq pour Tsurv ----------------
            ll=0
            if (Tsurv(i).eq.zi(nz)) then
               ll=nz-1
            end if
            som=0.d0
            do kk=2,nz
               if ((Tsurv(i).ge.zi(kk-1)).and.(Tsurv(i).lt.zi(kk))) &
                    then
                  ll=kk-1
               end if
            end do
            if (ll.gt.1) then
               do ii=1,ll-1
                  som=som+brisq(ii)
               end do
            end if

            surv(g)=som+brisq(ll)*Tim3(i)+brisq(ll+1)*Tim2(i) &
                 +brisq(ll+2)*Tim1(i)+brisq(ll+3)*Tim(i)
            risq(g)=brisq(ll)*Tmm3(i)+brisq(ll+1)*Tmm2(i)     &
                 +brisq(ll+2)*Tmm1(i)+brisq(ll+3)*Tmm(i)

            !------------ survie et risq pour Tsurv0 ----------------

            if (idtrunc.eq.1) then
               ll=0
               if (Tsurv0(i).eq.zi(nz)) then
                  ll=nz-1
               end if
               som=0.d0
               do kk=2,nz
                  if ((Tsurv0(i).ge.zi(kk-1)).and.(Tsurv0(i).lt.zi(kk))) &
                       then
                     ll=kk-1
                  end if
               end do
!               if (ll.lt.1.or.ll.gt.nz-1) then
!                  write(*,*) 'probleme dans fct_risq splines'
!                  write(*,*) 'll=',ll,'T=',Tsurv0(i)
!                  stop
!               end if
               if (ll.gt.1) then
                  do ii=1,ll-1
                     som=som+brisq(ii)
                  end do
               end if

               surv0(g)=som+brisq(ll)*Tim03(i)+brisq(ll+1)*Tim02(i) &
                    +brisq(ll+2)*Tim01(i)+brisq(ll+3)*Tim0(i)

            end if

            !------------ survie et risq pour Tsurvint ----------------


            if (ind_survint(i).eq.1) then

!            write(*,*)'i',i,tsurvint(i),ind_survint(i),tsurv(i),tsurv0(i)
!            write(*,*)timt3(i),Timt2(i),timt1(i)

               ll=0
               if (Tsurvint(i).eq.zi(nz)) then
                  ll=nz-1
               end if
               som=0.d0
               do kk=2,nz
                  if((Tsurvint(i).ge.zi(kk-1)).and.(Tsurvint(i).lt.zi(kk))) &
                       then
                     ll=kk-1
                  end if
               end do
!               if (ll.lt.1.or.ll.gt.nz-1) then
!                  write(*,*) 'probleme dans fct_risq splines'
!                  write(*,*) 'll=',ll,'T=',Tsurvint(i)
!                  stop
!               end if
               if (ll.gt.1) then
                  do ii=1,ll-1
                     som=som+brisq(ii)
                  end do
               end if

               survint(g)=som+brisq(ll)*Timt3(i)+brisq(ll+1)*Timt2(i) &
                    +brisq(ll+2)*Timt1(i)+brisq(ll+3)*Timt(i)

            else
               survint(g)=surv(g)
            end if

         end if


      end subroutine fct_risq_i



!---------------------------------------------------------------
!
      !SPLINES
!
!---------------------------------------------------------------







      subroutine splines()
     use commun_joint
    use commun_modele_joint
          implicit none

    integer::i,k,n,l
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
        hn,hh2,hh3



        l=0
        Tmm=0.d0
        Tmm1=0.d0
        Tmm2=0.d0
        Tmm3=0.d0
        Tim=0.d0
        Tim1=0.d0
        Tim2=0.d0
        Tim3=0.d0
        Tmm0=0.d0
        Tmm01=0.d0
        Tmm02=0.d0
        Tmm03=0.d0
        Tim0=0.d0
        Tim01=0.d0
        Tim02=0.d0
        Tim03=0.d0
        Tmmt=0.d0
        Tmmt1=0.d0
        Tmmt2=0.d0
        Tmmt3=0.d0
        Timt=0.d0
        Timt1=0.d0
        Timt2=0.d0
        Timt3=0.d0





        zi(-2)=zi(1)
    zi(-1)=zi(1)
    zi(0)=zi(1)
    zi(nz+1)=zi(nz)
        zi(nz+2)=zi(nz)
    zi(nz+3)=zi(nz)
        n=nz+2
!------------------- Tsurv ---------------------------
      Do i=1,ns

    do k=2,n-2
       if ((Tsurv(i).ge.zi(k-1)).and.  &
               Tsurv(i).lt.zi(k)) then
               l=k-1
           end if
    end do

    if (Tsurv(i).eq.zi(n-2)) then
       l=n-3
    end if

    ht = Tsurv(i)-zi(l)
    htm = Tsurv(i)-zi(l-1)
    h2t = Tsurv(i)-zi(l+2)
    ht2 = zi(l+1)-Tsurv(i)
    ht3 = zi(l+3)-Tsurv(i)
    hht = Tsurv(i)-zi(l-2)
    h = zi(l+1)-zi(l)
    hh = zi(l+1)-zi(l-1)
    h2 = zi(l+2)-zi(l)
    h3 = zi(l+3)-zi(l)
    h4 = zi(l+4)-zi(l)
    h3m = zi(l+3)-zi(l-1)
    h2n = zi(l+2)-zi(l-1)
    hn = zi(l+1)-zi(l-2)
    hh3 = zi(l+1)-zi(l-3)
    hh2 = zi(l+2)-zi(l-2)

    if (Tsurv(i).ne.zi(n-2)) then

       Tmm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
       Tmm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

    end if

    if (Tsurv(i).eq.zi(n-2)) then

       Tmm3(i) = 0.d0
       Tmm2(i) = 0.d0
       Tmm1(i) = 0.d0
       Tmm(i) = 4.d0/h

    end if

    Tim3(i) = (0.25d0*(Tsurv(i)-zi(l-3))*Tmm3(i)) &
                         +(0.25d0*hh2*Tmm2(i))        &
            +(0.25d0*h3m*Tmm1(i))+(0.25d0*h4*Tmm(i))
        Tim2(i) = (0.25d0*hht*Tmm2(i))  &
            +(h3m*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
        Tim1(i) = (htm*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
    Tim(i) = ht*Tmm(i)*0.25d0

!------------------- Tsurv0 --------------------------

        if (idtrunc.eq.1) then

           do k=2,n-2
              if ((Tsurv0(i).ge.zi(k-1)).and.   &
                   Tsurv0(i).lt.zi(k)) then
                 l=k-1
              end if
           end do

           if (Tsurv0(i).eq.zi(n-2)) then
              l=n-3
           end if

           ht = Tsurv0(i)-zi(l)
           htm = Tsurv0(i)-zi(l-1)
           h2t = Tsurv0(i)-zi(l+2)
           ht2 = zi(l+1)-Tsurv0(i)
           ht3 = zi(l+3)-Tsurv0(i)
           hht = Tsurv0(i)-zi(l-2)
           h = zi(l+1)-zi(l)
           hh = zi(l+1)-zi(l-1)
           h2 = zi(l+2)-zi(l)
           h3 = zi(l+3)-zi(l)
           h4 = zi(l+4)-zi(l)
           h3m = zi(l+3)-zi(l-1)
           h2n = zi(l+2)-zi(l-1)
           hn = zi(l+1)-zi(l-2)
           hh3 = zi(l+1)-zi(l-3)
           hh2 = zi(l+2)-zi(l-2)

           if (Tsurv0(i).ne.zi(nz-2)) then

              Tmm03(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

              Tmm02(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                   +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                   +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
              Tmm01(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                   +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                   +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
              Tmm0(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

           end if

           if (Tsurv0(i).eq.zi(n-2)) then

              Tmm03(i) = 0.d0
              Tmm02(i) = 0.d0
              Tmm01(i) = 0.d0
              Tmm0(i) = 4.d0/h

           end if

           Tim03(i) = (0.25d0*(Tsurv0(i)-zi(l-3))*Tmm03(i))  &
                +(0.25d0*hh2*Tmm02(i))           &
                +(0.25d0*h3m*Tmm01(i))+(0.25d0*h4*Tmm0(i))
           Tim02(i) = (0.25d0*hht*Tmm02(i))                  &
                +(h3m*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
           Tim01(i) = (htm*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
           Tim0(i) = ht*Tmm0(i)*0.25d0

        end if


!------------------- Tsurvint --------------------------
        if (ind_survint(i).eq.1) then
        do k=2,n-2
          if ((Tsurvint(i).ge.zi(k-1)).and. &
                 Tsurvint(i).lt.zi(k)) then
                  l=k-1
              end if
       end do

        if (Tsurvint(i).eq.zi(nz-2)) then
          l=n-3
       end if

       ht = Tsurvint(i)-zi(l)
       htm = Tsurvint(i)-zi(l-1)
       h2t = Tsurvint(i)-zi(l+2)
       ht2 = zi(l+1)-Tsurvint(i)
       ht3 = zi(l+3)-Tsurvint(i)
       hht = Tsurvint(i)-zi(l-2)
       h = zi(l+1)-zi(l)
       hh = zi(l+1)-zi(l-1)
       h2 = zi(l+2)-zi(l)
       h3 = zi(l+3)-zi(l)
       h4 = zi(l+4)-zi(l)
       h3m = zi(l+3)-zi(l-1)
       h2n = zi(l+2)-zi(l-1)
       hn = zi(l+1)-zi(l-2)
       hh3 = zi(l+1)-zi(l-3)
       hh2 = zi(l+2)-zi(l-2)

       if (Tsurvint(i).ne.zi(nz-2)) then

          Tmmt3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
          Tmmt2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn)) &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))     &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
              Tmmt1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))       &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
              Tmmt(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

       end if

       if (Tsurvint(i).eq.zi(nz-2)) then

          Tmmt3(i) = 0.d0
          Tmmt2(i) = 0.d0
          Tmmt1(i) = 0.d0
          Tmmt(i) = 4.d0/h

       end if

       Timt3(i) = (0.25d0*(Tsurvint(i)-zi(l-3))*Tmmt3(i)) &
                         +(0.25d0*hh2*Tmmt2(i))               &
            +(0.25d0*h3m*Tmmt1(i))+(0.25d0*h4*Tmmt(i))
           Timt2(i) = (0.25d0*hht*Tmmt2(i))                   &
             +(h3m*Tmmt1(i)*0.25d0)+(h4*Tmmt(i)*0.25d0)
           Timt1(i) = (htm*Tmmt1(i)*0.25d0)+(h4*Tmmt(i)*0.25d0)
       Timt(i) = ht*Tmmt(i)*0.25d0

        else
           Timt3(i) =Tim3(i)
           Timt2(i) =Tim2(i)
           Timt1(i) =Tim1(i)
           Timt(i) =Tim(i)
        end if

       End Do

      end subroutine splines


!-------------------------------------------------------------------------
!
      !SCORE TEST
!
!-------------------------------------------------------------------------





      subroutine scoretest(b,npm,statsc)
     use commun_joint
    use commun_modele_joint
      use optim
      implicit none

      double precision,intent(out)::statsc

      double precision::survmarg,risqcum,risqcum0,coef
      double precision,dimension(nea,nea)::totretard,retard,varSi
      double precision,dimension(:,:),allocatable::stat,ppi
      double precision,dimension(nea):: statscore,RE
!      double precision,dimension(ns):: Somtest
      double precision,dimension(nea,nea)::varscore,varscore2
      double precision,dimension(nea*(nea+1)/2)::vars,vars2
      double precision,dimension(nea,nea)::BB,Ebb,Varb
      double precision,dimension(nea,maxmes)::BZT,BZV
      double precision,dimension(maxmes,nea)::  Z,P
      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk &
      ,nef,iqq
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Se
      double precision,dimension(npm) :: b,b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision ::Y4,f,fevt
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3

      double precision,dimension(nvarxevt)::bevt,Xevt
      double precision::bevtint
      double precision,dimension(nprisq)::brisq
!     double precision,dimension(ns,ng)::risq,surv,surv0, &
!      survint,risq2
      double precision,dimension(:,:),allocatable::risq,surv,surv0, &
      survint,risq2

      eps=1.d-20
      b1=0.d0
      do k=1,npm
         b1(k)=b(k)
      end do



      allocate(risq(ns,ng),surv(ns,ng),surv0(ns,ng),survint(ns,ng),risq2(ns,ng))
      allocate(stat(ns,nea),ppi(ns,ng))
      stat=0.d0
      ppi=0.d0
      risq=0.d0
      surv=0.d0
      surv0=0.d0
      survint=0.d0

      allocate(VC(maxmes,maxmes),Se(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      Vi=0.d0
      Se=0.d0
      risq2=0.d0
      statsc=0.d0



!----------- rappel des parametres utilises ---------

      nef=nprob+nrisq+nvarxevt+ncssg+ncg*ng

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

!------------------- calcul de risq ------------------------
! calcul de brisq e chaque composante et risq, surv et surv0 pour tous les i et tous les g


      risq=0.d0
      surv=0.d0
      surv0=0.d0
      survint=0.d0

      do g=1,ng

         brisq=0.d0
         if (logspecif.eq.1.or.typrisq.eq.2) then
               if (risqcom.eq.0) then
                 do k=1,nprisq
                   brisq(k)=exp(b1(nprob+nprisq*(g-1)+k))
               end do
               elseif (risqcom.eq.1) then
                do k=1,nprisq
                  brisq(k)=exp(b1(nprob+k))
                end do
               elseif (risqcom.eq.2) then
             do k=1,nprisq
                brisq(k)=exp(b1(nprob+k))
              end do
            endif
        else
               if (risqcom.eq.0) then
              do k=1,nprisq
                  brisq(k)=b1(nprob+nprisq*(g-1)+k)   &
                      *b1(nprob+nprisq*(g-1)+k)
               end do
               elseif (risqcom.eq.1) then
                do k=1,nprisq
                  brisq(k)=b1(nprob+k) *b1(nprob+k)

                end do
               elseif (risqcom.eq.2) then
               do k=1,nprisq
                brisq(k)=b1(nprob+k)*b1(nprob+k)
               end do
              end if
        end if
!         write(*,*)'g=',g,'brisq=',(brisq(k),k=1,nrisq)

         call fct_risq(brisq,g,risq,surv,surv0,survint)

         if (risqcom.eq.2.and.ng.gt.1.and.g.lt.ng) then
            do i=1,ns
               risq(i,g)=risq(i,g)*exp(b1(nprob+nprisq+g))
               surv(i,g)=surv(i,g)*exp(b1(nprob+nprisq+g))
               survint(i,g)=survint(i,g)*exp(b1(nprob+nprisq+g))
               surv0(i,g)=surv0(i,g)*exp(b1(nprob+nprisq+g))
            end do
         end if

      end do


      varscore=0.d0
      statscore=0.d0
      totretard=0.d0
      stat=0.d0

! ----------- boucle sur les individus -------------
      kk=0
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


         Xevt=0.d0
            l=0
            do k=1,nv
               if (idxevt(k).eq.1.or.idxevt(k).eq.2) then
                  l=l+1
                  Xevt(l)=X(it+1,k)
               end if

            end do
! creation de s2*I et Y1

         Se=0.d0
     Y1=0.d0
         do j=1,nmes(i)
        kk=kk+1
            Se(j,j)=b1(npm)*b1(npm)
        Y1(j)=dble(Y(kk))
         end do

! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

        if (nwg.eq.0) then

          P=0.d0
      P=MATMUL(Z,Ut)
      BZT=0.d0
      BZT=MATMUL(Ut,transpose(P))
      BB=0.d0
      BB=MATMUL(Ut,transpose(Ut))
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

         CALL DSINV(Vi,nmes(i),eps,ier,det)
         if (ier.eq.-1) then
            statsc=9999.d0
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



! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       if (prior(i).ne.0) then
           pi=0.d0
           pi(prior(i))=1.d0
       else
           Xprob=0.d0
           Xprob(1)=1
           l=0
           do k=1,nv
              if (idprob(k).eq.1) then
                 l=l+1
                 Xprob(1+l)=X(it+1,k)
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
       retard=0.d0
       survmarg=0.d0
       varSi=0.d0
       b0=0.d0
       b2=0.d0

       do g=1,ng
          nmoins=0
          l2=0
          m2=0
          do k=1,nv
             if (idg(k).eq.1) then
                m2=m2+1
                b0(m2)=b1(nprob+nrisq+nvarxevt+nmoins+1)
                nmoins=nmoins+1
             else if (idg(k).eq.2) then
                l2=l2+1
                b2(l2)=b1(nprob+nrisq+nvarxevt+nmoins+g)
                nmoins=nmoins+ng
             end if
          end do

! bevt
             bevt=0.d0
            bevtint=0.d0
               m=0
               l=1
               do k=1,nv
                  if (idxevt(k).eq.1) then
                     bevt(l)=b1(nprob+nrisq+m+1)
                     l=l+1
                     m=m+1
                  end if
                  if (idxevt(k).eq.2) then
                     bevt(l)=b1(nprob+nrisq+m+g)
                     l=l+1
                     m=m+ng
                  end if

               end do
               if (nvdepsurv.ne.0) then
                  bevtint=b1(nprob+nrisq+nvarxevt)
               end if

!    variance covariance si spec aux classes :

          if (nwg.ne.0) then
             Ut1=0.d0
             if (g.eq.ng) then
                Ut1=Ut
             else
                Ut1=Ut*abs(b1(nef+nvc+g))
             end if

             P=0.d0
         P=MATMUL(Z,Ut1)
         BZT=0.d0
         BZT=MATMUL(Ut1,transpose(P))
         BB=0.d0
         BB=MATMUL(Ut1,transpose(Ut1))
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

             CALL DSINV(Vi,nmes(i),eps,ier,det)
             if (ier.eq.-1) then
                statsc=9999.d0
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


!- ------------------------------
          mu=0.d0
      mu=matmul(X0,b0)+matmul(X2,b2)
      Y2=Y1-mu
      Y3=MATMUL(VC,Y2)
      Y4=DOT_PRODUCT(Y2,Y3)
      RE=0.d0
      RE=MATMUL(BZT,Y3)
      BZV=MATMUL(BZT,VC)
      varb=0.d0
      varb=BB-MATMUL(BZV,transpose(BZT))
      Ebb=0.d0

          do j=1,nea
             do k=1,nea
              Ebb(j,k)=Varb(j,k)+RE(j)*RE(k)
             end do
          end do

!          write(*,*)i,'RE, varb',RE(1),varb(1,1)

      risqcum=0.d0
      risqcum0=0.d0
!--------- calcul de fevt --------------
          risq2(i,g)=risq(i,g)
          if (ind_survint(i).eq.1) then
             risq2(i,g)=risq(i,g)*exp(bevtint)
          end if
          fevt=log(risq2(i,g))+DOT_PRODUCT(bevt,Xevt)

          fi(g)=fi(g)- nmes(i)*log(dble(2*3.14159265))
          fi(g)=fi(g)-det
          fi(g)=fi(g)- Y4
          fi(g)=fi(g)/(2.d0)

          risqcum=exp(DOT_PRODUCT(Xevt,bevt))*               &
                 (survint(i,g)+dexp(bevtint)*(surv(i,g)   &
                 -survint(i,g)))
!       write(*,*)'i',i,fi(g),Devt(i),fevt,risqcum 
           temp=fi(g)+Devt(i)*fevt-risqcum
       fi(g)=0.d0     
       fi(g)=dexp(temp)

! statr du score et var pour sujet i :test global


!          write(*,*)i,'ebb',((Ebb(j,k),j=1,nea),k=1,nea)
         do j=1,nea
             stat(i,j)=stat(i,j)+pi(g)*fi(g)*(Devt(i)-risqcum)*RE(j)
             do k=1,nea
                varSi(j,k)=varSi(j,k)+pi(g)*fi(g)*  &
               ((Devt(i)-risqcum)**2-risqcum)*Ebb(j,k)
             end do
          end do

! correction variance si troncature et si test intercept + pente

         if(idtrunc.eq.1) then
            risqcum0=exp(DOT_PRODUCT(Xevt,bevt))*surv0(i,g)
            survmarg=survmarg+pi(g)*exp(-risqcum0)
            coef=pi(g)*exp(-risqcum0) &
                             *risqcum0*(risqcum0-1.d0)
            retard=retard+coef*BB
         end if

!          write(*,*)'i',i,g,risqcum,risqcum0
       end do !Fin boucle sur les groupes g

! contribution vrais du  sujet i=f(Yi,Ti,Di))

       f=DOT_PRODUCT(pi,fi)

! stat et var du score
       do j=1,nea
          stat(i,j)=stat(i,j)/f
          statscore(j)= statscore(j)+stat(i,j)
       end do

!       write(*,*)i,'VarSi',((varSi(j,k),j=1,nea),k=1,nea)
       do j=1,nea
          do k=1,nea
             varscore(j,k)=varscore(j,k)-varSi(j,k)/f + (stat(i,j)*stat(i,k))
             if (idtrunc.eq.1) then
                 totretard(j,k)=totretard(j,k) + retard(j,k)/survmarg
             end if
          end do
       end do
! verification de la somme des probabilites egalent a 1:Somtest
!       Somtest=0.d0
       do g=1,ng
      ppi(i,g)=pi(g)*fi(g)/f
!      Somtest(i)=Somtest(i)+ppi(i,g)
       end do

    it=it+nmes(i)
    end do
! fin boucle sujets

! Correction variance si troncature
      if(idtrunc.eq.1) then
         varscore=varscore+totretard
      end if


! Calcul de la variance empirique du score (cf Freedman Am Stat 2007)

       varscore2=0.D0
       DO j=1,nea
        DO k=1,nea
          Do i=1,ns
            varscore2(j,k)=varscore2(j,k)+ Stat(i,j)*  &
     Stat(i,k)
          END DO
          varscore2(j,k)=varscore2(j,k)- Statscore(j)  &
     *Statscore(k)/dble(ns)
        END DO
       END DO


       vars=0.d0
       vars2=0.d0

      jj=0
      do j=1,nea
         do k=j,nea
            jj=j+k*(k-1)/2
            Vars(jj)=Varscore(j,k)
            Vars2(jj)=Varscore2(j,k)
         end do
      end do

!      CALL DSINV(vars,nea,eps,ier,det)
!      if (ier.eq.-1) then
!         write(*,*)'dsinv varscore ier',ier,'det',det
!         write(*,*)'VarsCore', (Vars(j),j=1,nea*(nea+1)/2)
!      endif

      CALL DSINV(vars2,nea,eps,ier,det)
      if (ier.eq.-1) then
!         write(*,*)'dsinv varscore2 ier',ier,'det',det
!         write(*,*)'VarsCore2', (Vars2(j),j=1,nea*(nea+1)/2)
         statsc=9999.d0
         goto 147
      else
         vars=vars2
      end if

      statsc=0.D0

      DO j=1,nea
         DO k=1,nea
            if(k.ge.j) then
             iqq=k*(k-1)/2+j
            else
             iqq=j*(j-1)/2+k
            end if
           statsc=statsc+statscore(j)*vars(iqq)*statscore(k)
         END DO
      END DO


 147  continue
      deallocate(VC,Se,Vi,stat,ppi)
      deallocate(risq,surv,surv0,survint,risq2)


      if (statsc.ne.statsc) statsc=9999.d0

!    write(*,*)'statsc',statsc
      return

      end subroutine scoretest





!-------------------------------------------------------------------
! IS NAN
!-------------------------------------------------------------------

      logical function is_nan(x)

      implicit none

      double precision,intent(in)::x

      if (x .ne. x) then
       is_nan=.true.
      else
           is_nan=.false.
      end if

      return

       end




!-------------------------------------------------------------------
!
!!  FCT_RISQ ESTIME
!
!-------------------------------------------------------------------




      subroutine fct_risq_estime(brisq,time,nsim,g,risq,surv)
    use commun_joint
    use commun_modele_joint
 !     use communJ

      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


      double precision, dimension(nprisq)::brisq
      integer ::i,g,ll,kk,ii,l,j,nsim
      double precision,dimension(nsim*ng):: risq,surv
      double precision::som
      double precision,dimension(nsim)::time


      l=0
      do i=1,nsim
         if (typrisq.eq.2) then
            surv(nsim*(g-1)+i)=brisq(1)*(time(i))**brisq(2)
            risq(nsim*(g-1)+i)=brisq(1)*brisq(2) &
                 *(time(i))**(brisq(2)-1)
         end if

         if (typrisq.eq.1) then
            do j=1,nz-1
               som=0.d0
               do l=1,j-1
                  som=som+brisq(l)*(zi(l+1)-zi(l))
               end do
               if (time(i).ge.zi(j).and.time(i).le.zi(j+1)) then
                  surv(nsim*(g-1)+i)=som+brisq(j)*(time(i)-zi(j))
                  risq(nsim*(g-1)+i)=brisq(j)
               end if
            end do
         end if

         if (typrisq.eq.3) then
            !------------ survie et risq pour Tsurv ----------------
            ll=0
            if (time(i).eq.zi(nz)) then
               ll=nz-1
            end if
            som=0.d0
            do kk=2,nz
               if ((time(i).ge.zi(kk-1)).and.(time(i).lt.zi(kk))) &
                    then
                  ll=kk-1
               end if
            end do
            if (ll.gt.1) then
               do ii=1,ll-1
                  som=som+brisq(ii)
               end do
            end if

            surv(nsim*(g-1)+i)=som+brisq(ll)*Tim3_est(i)+brisq(ll+1)*Tim2_est(i) &
                 +brisq(ll+2)*Tim1_est(i)+brisq(ll+3)*Tim_est(i)
            risq(nsim*(g-1)+i)=brisq(ll)*Tmm3_est(i)+brisq(ll+1)*Tmm2_est(i)     &
                 +brisq(ll+2)*Tmm1_est(i)+brisq(ll+3)*Tmm_est(i)

         end if
      end do



    end subroutine fct_risq_estime







!-------------------------------------------------------------------
!
!!  SPLINES ESTIMES
!
!-------------------------------------------------------------------






      subroutine splines_estime(time,nsim)
    use commun_joint
    use commun_modele_joint
          implicit none

    integer::i,k,n,l
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
        hn,hh2,hh3
        double precision,dimension(nsim)::time
        integer::nsim


        Tmm_est=0.d0
        Tmm1_est=0.d0
        Tmm2_est=0.d0
        Tmm3_est=0.d0
        Tim_est=0.d0
        Tim1_est=0.d0
        Tim2_est=0.d0
        Tim3_est=0.d0

        l=0


        zi(-2)=zi(1)
        zi(-1)=zi(1)
        zi(0)=zi(1)
        zi(nz+1)=zi(nz)
        zi(nz+2)=zi(nz)
        zi(nz+3)=zi(nz)
        n=nz+2
!------------------- time ---------------------------
      Do i=1,nsim
          do k=2,n-2
       if ((time(i).ge.zi(k-1)).and.  &
               time(i).lt.zi(k)) then
               l=k-1
           end if
    end do

    if (time(i).eq.zi(n-2)) then
       l=n-3
    end if

    ht = time(i)-zi(l)
    htm = time(i)-zi(l-1)
    h2t = time(i)-zi(l+2)
    ht2 = zi(l+1)-time(i)
    ht3 = zi(l+3)-time(i)
    hht = time(i)-zi(l-2)
    h = zi(l+1)-zi(l)
    hh = zi(l+1)-zi(l-1)
    h2 = zi(l+2)-zi(l)
    h3 = zi(l+3)-zi(l)
    h4 = zi(l+4)-zi(l)
    h3m = zi(l+3)-zi(l-1)
    h2n = zi(l+2)-zi(l-1)
    hn = zi(l+1)-zi(l-2)
    hh3 = zi(l+1)-zi(l-3)
    hh2 = zi(l+2)-zi(l-2)

    if (time(i).ne.zi(n-2)) then

       Tmm3_est(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
       Tmm2_est(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm1_est(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm_est(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

    end if

    if (time(i).eq.zi(n-2)) then

       Tmm3_est(i) = 0.d0
       Tmm2_est(i) = 0.d0
       Tmm1_est(i) = 0.d0
       Tmm_est(i) = 4.d0/h

    end if
    
    Tim3_est(i) = (0.25d0*(time(i)-zi(l-3))*Tmm3_est(i)) &
         +(0.25d0*hh2*Tmm2_est(i))        &
         +(0.25d0*h3m*Tmm1_est(i))+(0.25d0*h4*Tmm_est(i))
    Tim2_est(i) = (0.25d0*hht*Tmm2_est(i))  &
         +(h3m*Tmm1_est(i)*0.25d0)+(h4*Tmm_est(i)*0.25d0)
    Tim1_est(i) = (htm*Tmm1_est(i)*0.25d0)+(h4*Tmm_est(i)*0.25d0)
    Tim_est(i) = ht*Tmm_est(i)*0.25d0
    
 End Do
 
end subroutine splines_estime
