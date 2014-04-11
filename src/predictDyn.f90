
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   Calcul de la proba de survenue d'un evenement dans une fenetre (s,s+t)  !
   !      a partir d'un modele conjoint                                        !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 subroutine proba_evt_dyn(Y0,X0,ns0,nmes0,nobs0,prior0,ng0,nv0,idiag0,nwg0,ncor0,logspecif0, &
    idxevt0,zi0,idea0,idg0,idprob0,idcor0,risqcom0,nvdepsurv0,nvarxevt0,typrisq0,nz0, &
    tsurv00,tsurvint0,idtrunc0,b,npm0,timemes,landmark,horizon,nbland,nbhoriz,proba)
    
    use commun_joint
    use commun_modele_joint

    IMPLICIT NONE


    !Declaration des variables en entree
    integer,intent(in)::nobs0,ns0,nv0,ng0,nz0,npm0,nbland,nbhoriz
    double precision,dimension(nobs0),intent(in):: Y0
    double precision,dimension(nobs0*nv0),intent(in)::X0
    integer,dimension(ns0),intent(in)::nmes0,prior0
    double precision,dimension(ns0),intent(in)::tsurv00,tsurvint0
    integer,intent(in)::idiag0,nwg0,ncor0,logspecif0,risqcom0,nvdepsurv0,nvarxevt0,typrisq0,idtrunc0
    integer,dimension(nv0),intent(in)::idxevt0,idea0,idg0,idprob0,idcor0
    double precision,dimension(nz0),intent(in)::zi0
    double precision,dimension(npm0),intent(in)::b
    double precision,dimension(nobs),intent(in)::timemes
    double precision,dimension(nbland),intent(in)::landmark
    double precision,dimension(nbhoriz),intent(in)::horizon
    
    ! Sortie
    double precision,dimension(ns0*nbland*nbhoriz),intent(out)::proba
    
    ! Variables locales
    integer::k,g,i,l,m,kk,sumnmes,j,ks,kt,ist    !,nsim
    double precision,dimension(ns0,ng0)::ppi
    double precision,dimension(nz0+2)::brisq   !dim surestimee
    !double precision,dimension(2*ng0)::risq,surv
    double precision,dimension(ns0,ng0)::risq,surv1,surv2,surv0,survint1,survint2
    double precision,dimension(nvarxevt0)::bevt,Xevt  ! surestime
    double precision,dimension(1)::bevtint
    
    ! Allocation de memoire
    allocate(Y(nobs0),X(nobs0,nv0),nmes(ns0),idxevt(nv0),idea(nv0), &
    idg(nv0),idprob(nv0),idcor(nv0),prior(ns0),tsurv(ns0),tsurv0(ns0), &
    tsurvint(ns0),ind_survint(ns0))
    
    select case (typrisq0)
        case (1)
          allocate(zi(nz0))
          zi(1:nz0)=zi0(1:nz0)
        case (2)
          allocate(zi(nz0))
          zi(1:nz0)=zi0(1:nz0)
        case (3)
           allocate(zi(-2:nz0+3))
           do i=1,nz0
              zi(i)=zi0(i) ! reste rempli dans subroutine splines
           end do
    end select
                          
      ! print*,"nmes0",nmes0,"ns0=",ns0
    ! Remplir les modules avec les arguments
    nobs=nobs0
    ns=ns0
    nv=nv0
    ng=ng0
    nz=nz0
    Y(1:nobs)=Y0(1:nobs)
    prior(1:ns)=prior0(1:ns)
    idiag=idiag0
    nwg=nwg0
    ncor=ncor0
    logspecif=logspecif0
    risqcom=risqcom0
    nvdepsurv=nvdepsurv0
    nvarxevt=nvarxevt0
    tsurv0(1:ns)=tsurv00(1:ns)
    idtrunc=idtrunc0
    typrisq=typrisq0
    idxevt(1:nv)=idxevt0(1:nv)
    idea(1:nv)=idea0(1:nv)
    idg(1:nv)=idg0(1:nv)
    idprob(1:nv)=idprob0(1:nv)
    idcor(1:nv)=idcor0(1:nv)
    do k=1,nv
     X(1:nobs,k)=X0(((k-1)*nobs+1):(k*nobs))
    end do
          
    ! Remplir les autres variables des modules
   
    nprisq=0
    if (typrisq.eq.1) then
      nprisq=nz-1
    end if
    if (typrisq.eq.2) then
      nprisq=2
    end if
    if (typrisq.eq.3) then
      nprisq=nz+2
    end if
        
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
    
   
    nea=0
    ncg=0
    ncssg=0
    nprob=ng-1  ! car indicateur pour intercept dans idprob vaut tjrs 0 mais intercept tjrs present
    nvarprob=min(ng-1,1)
    do k=1,nv
       if (idg(k).eq.1) then
          ncssg=ncssg+1      ! nb var. sans melange
       else if (idg(k).eq.2) then
          ncg=ncg+1          ! nb var. dans melange
       end if
       nea=nea+idea(k)
       nprob=nprob+(idprob(k))*(ng-1)
       nvarprob=nvarprob+idprob(k)
    end do

    !  nb effets fixes = nb effets fixes sans melange
    !                  + ng fois le nb de var dans melange

    nvc=0
    if (idiag.eq.1) then
       nvc=nea
    else if(idiag.eq.0) then
       nvc=nea*(nea+1)/2
    end if
    
    maxmes=maxval(nmes0)
    
    if (typrisq.eq.3) then
     
      allocate(Tmm(ns0),Tmm1(ns0),Tmm2(ns0),Tmm3(ns0),Tim(ns0)         &
     ,Tim1(ns0),Tim2(ns0),Tim3(ns0),Tmm0(ns0),Tmm01(ns0),Tmm02(ns0)  &
     ,Tmm03(ns0),Tim0(ns0),Tim01(ns0),Tim02(ns0),Tim03(ns0),          &
      Tmmt(ns0),Tmmt1(ns0),Tmmt2(ns0),Tmmt3(ns0),Timt(ns0),Timt1(ns0) &
     ,Timt2(ns0),Timt3(ns0))
      
    end if
    
    ! modules ok (sauf tsurv,nmes,Tmm,tsurvint,ind_survint  a faire pour chaque landmark)

    proba=0.d0   
    
    do ks=1,nbland
      ! print*,"ks=",ks
      tsurv(1:ns0)=landmark(ks)
      tsurvint(1:ns0)=landmark(ks)
      ind_survint(1:ns0)=0
      
      nmes=0
      sumnmes=0
      do i=1,ns
        do j=1,nmes0(i)
        !print*,"timemes",timemes(sumnmes+j),i,j
          if(timemes(sumnmes+j)<landmark(ks)) then
            nmes(i) = nmes(i) + 1 
          end if
        end do
        !   print*,"nmes=",nmes

        
        if(tsurvint0(i)<landmark(ks)) then
          tsurvint(i) = tsurvint0(i)
          ind_survint(i) = 1
        end if      
        
        sumnmes = sumnmes+nmes0(i)
      end do  
    
      if (typrisq.eq.3) then
        call splines() ! pour remplir Tmm etc
      end if
       
      ! calcul des proba a posteriori
    
      ppi=0.d0
    
      !if (ng.gt.1) then
        call postprobj2(b,npm0,nmes0,ppi)      
 !         print*,ppi
      !else
      !  ppi=1.d0
      !end if

      ! ok proba dans ppi
    
    
      risq=0.d0
      surv1=0.d0
      surv0=0.d0
      survint1=0.d0    
      surv2=0.d0
      survint2=0.d0
               
      do g=1,ng
      
      ! calcul du risque cumule (Lambda0) aux temps s 
      
        brisq=0.d0
       
        if (logspecif.eq.1) then
        
          if (risqcom.eq.0) then
            do k=1,nprisq
              brisq(k)=exp(b(nprob+nprisq*(g-1)+k))
            end do
          elseif (risqcom.eq.1) then
            do k=1,nprisq
              brisq(k)=exp(b(nprob+k))
            end do
          elseif (risqcom.eq.2) then
            do k=1,nprisq
              brisq(k)=exp(b(nprob+k))    ! contient le coef de proportionnalite
            end do
          end if
          
        else   ! ie  pas de logspecif 
        
          if (risqcom.eq.0) then
            do k=1,nprisq
              brisq(k)=b(nprob+nprisq*(g-1)+k)*b(nprob+nprisq*(g-1)+k)
            end do
          elseif (risqcom.eq.1) then
            do k=1,nprisq
              brisq(k)=b(nprob+k)*b(nprob+k)
            end do
          elseif (risqcom.eq.2) then
            do k=1,nprisq
              brisq(k)=b(nprob+k)*b(nprob+k)
            end do
          end if
          
        end if
 !      print*,brisq
        tsurv=landmark(ks) 
        
        do i=1,ns
          if(ind_survint(i).eq.0) then
            tsurvint(i) = landmark(ks)
          end if
        end do  
        
        if (typrisq.eq.3) then  
          call splines() ! actualiser Tmm etc
        end if        
        
        call fct_risq(brisq,g,risq,surv1,surv0,survint1) 
         ! on a Lambda_0(s) de la classe g pour tous les sujets
!print*,'surv(s)',surv1
         
        ! print*,"tsurv=",tsurv
        ! print*,"Lambda(s)=",surv1
        ! print*,"tsurvint=",tsurvint
        ! print*,"Lambda(tt)",tsurvint

        ! creation de bevt (de la classe g) et bevtint
        bevt=0.d0
        bevtint=0.d0
        m=0
        l=0
        do k=1,nv
        
          if (idxevt(k).eq.1) then
             l=l+1        
             bevt(l)=b(nprob+nrisq+m+1)
             m=m+1
          end if
  
          if (idxevt(k).eq.2) then
             l=l+1        
             bevt(l)=b(nprob+nrisq+m+g)      
             m=m+ng
          end if
        end do
        
        if (nvdepsurv.ne.0) then
          bevtint(1)=b(nprob+nrisq+nvarxevt)
        end if
        ! bevt et bevtint ok
        
        do kt=1,nbhoriz
        
        ! calcul du risque cumule (Lambda0) aux temps s+t 
        
          tsurv=landmark(ks)+horizon(kt) 
          do i=1,ns
            if(ind_survint(i).eq.0) then
              tsurvint(i) = landmark(ks)+horizon(kt)
            end if
          end do  

          
          if (typrisq.eq.3) then  
            call splines() ! actualiser Tmm etc
          end if

          
          call fct_risq(brisq,g,risq,surv2,surv0,survint2)
            ! on a Lambda_0(s+t) de la classe g  pour tous les sujets
          
        ! print*,"tsurv=",tsurv
        ! print*,"Lambda(s+t)=",surv2
        ! print*,"tsurvint=",tsurvint
        ! print*,"Lambda(tt)",survint2  

          sumnmes=0
    
          do i=1,ns
          

            
            ! creation de Xevt
            Xevt=0.d0
            kk=0
            
            do k=1,nv
          
              if (idxevt(k).eq.1.or.idxevt(k).eq.2) then
                kk=kk+1
                Xevt(kk)=X(sumnmes+1,k) ! +1 ok car variables ne dependent pas du temps
              end if
            end do
            
            ! proba = sum_g ( ppi_g * (S_ig(s)) - S_ig(s+t))/S_ig(s)) 
            ! avec S_ig(s) = exp(-Lambda0_g(s)* exp(Xevt*bevt))
        
            ist = (i-1)*nbland*nbhoriz + (ks-1)*nbhoriz + kt ! indice de proba a remplir
            
       !     if(i.eq.1.and.ks.eq.1.and.kt.eq.2) then 
          !   if(ng.eq.1) then
          ! print*,"ietg:",i,g
       !     print*,"ist=",ist
           !print*,"ppi=",ppi(i,g)
       !     print*,"surv1=",surv1(i,g)
          !  print*,"survint1=",survint1(i,g)
       !     print*,"surv2=",surv2(i,g)
          !  print*,"survint2=",survint2(i,g)
          ! print*,"i=",i,"Xb=",exp(dot_product(Xevt,bevt))
          !      end if

!print*,"i=",i,"nmes0=",nmes0(i),"nmes=",nmes(i)
!                 print*,i,g,ppi(i,g),ist
            if(ppi(i,g)>-1) then       
                  
              proba(ist) = proba(ist) + ppi(i,g) * &
                 ( exp(-(survint1(i,g)+exp(bevtint(1))*(surv1(i,g)-survint1(i,g)))*exp(dot_product(Xevt,bevt))) - &
                   exp(-(survint2(i,g)+exp(bevtint(1))*(surv2(i,g)-survint2(i,g)))*exp(dot_product(Xevt,bevt))) ) / &
                 exp(-(survint1(i,g)+exp(bevtint(1))*(surv1(i,g)-survint1(i,g)))*exp(dot_product(Xevt,bevt)))             
            else
            
              proba(ist) = -5.d0
              
            end if  
              
            sumnmes = sumnmes + nmes0(i)
    
          end do ! fin boucle sur i
          
        end do  ! fin boucle sur kt
      
      end do   ! fin boucle sur g
      
    end do     ! fin boucle sur ks

   ! ok
           ! print*,"proba=",proba
           
    if(typrisq.eq.3) then
 !     deallocate(Tmm_est,Tmm1_est,Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est)
    
      deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim  &       
     ,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02    &
     ,Tmm03,Tim0,Tim01,Tim02,Tim03,      &    
      Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1  &
     ,Timt2,Timt3)
    end if
    
    deallocate(Y,X,nmes,idxevt,zi,idea,idg,idprob,idcor,prior,tsurv,tsurv0,tsurvint,ind_survint)
    
 end subroutine proba_evt_dyn
 
 
 
!-------------------------------------------------------------
!
!          Subroutine pour calculer les
!      probabilites a posteriori de suivre chacune
!      des composantes g pour chacun des sujets i
!      pour des donn\E9es longitudinales tronqu\E9es en s
!-------------------------------------------------------------
! -> si aucune mesure a des temps < s , alors ppi=-5

      subroutine postprobj2(b,npm,nmes0,ppi)

      use commun_joint
      use commun_modele_joint
      use optim ! pour dsinv et dmfsd

      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk,nef,j1,j2
      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Corr
      double precision,dimension(npm) :: b,b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision ::Y4
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(ns,ng) ::ppi
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      integer,dimension(ns)::nmes0

      double precision,dimension(nvarxevt)::bevt,Xevt
      double precision,dimension(ns,ng)::risq2
      double precision::bevtint
      double precision,dimension(nprisq)::brisq
      double precision,dimension(ns,ng)::risq,surv,surv0,  &
      survint

      eps=1.D-20
      ppi=-5.D0
      b1=0.d0
      do k=1,npm
         b1(k)=b(k)
      end do

      allocate(VC(maxmes,maxmes),Corr(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      Vi=0.d0
      Corr=0.d0
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
         if (logspecif.eq.1) then
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
               survint(i,g)=survint(i,g)*exp(b1(nprob+nprisq+g))  ! si logscale=TRUE ou tjrs????
               surv0(i,g)=surv0(i,g)*exp(b1(nprob+nprisq+g))
            end do
         end if

      end do

! ----------- boucle sur les individus -------------
      kk=0
      it=0
      do i=1,ns

       if(nmes(i)>0) then

! -------- creation de Vi = ZiGZi+Ri+se*seIni ----------
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
                        
!matrice Corri=Ri+s2*I
        
        Corr=0.d0
        tcor=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,nmes(i)
                    tcor(j) = X(it+j,k)
                 end do
              end if
           end do
         end if
         do j1=1,nmes(i)
            do j2=1,nmes(i)
               if (j1.eq.j2) then
                 Corr(j1,j2) = b1(npm)*b1(npm)
               end if
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*exp(-b1(npm-2)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do

                        
! creation de Y1

         Y1=0.d0
         do j=1,nmes(i)
            Y1(j)=dble(Y(it+j))
         end do

! creation de P=Zi*Ut et V=P*P que si non spec aux classes

        if (nwg.eq.0) then
          P=0.d0
          P=MATMUL(Z,Ut)
          VC=0.d0
          VC=MATMUL(P,transpose(P))+Corr

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
             VC=MATMUL(P,transpose(P))+Corr



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


          fi(g)=fi(g) - nmes(i)*log(dble(2*3.14159265))
          fi(g)=fi(g) - det
          fi(g)=fi(g) - Y4
          fi(g)=fi(g)/(2.d0)

          ! NB : on suppose que Devt=0 pour tout le monde

          fi(g)=fi(g)-exp(DOT_PRODUCT(Xevt,bevt))* &
              (survint(i,g)+exp(bevtint)*(surv(i,g)-survint(i,g)))

          fi(g)=exp(fi(g))

       end do

       do g=1,ng
          ppi(i,g)=pi(g)*fi(g)/DOT_PRODUCT(pi,fi)
       end do

      end if ! fin if(nmes(i)>0)
       
       it=it+nmes0(i)

      end do


 147  continue


      deallocate(VC,Vi,Corr)

      return

      end subroutine postprobj2 
    
