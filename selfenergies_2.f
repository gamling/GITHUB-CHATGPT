************************************************************************
*****************************SIGMAF_HEAVY_U*****************************
************************************************************************
      SUBROUTINE SGMAFU(beta,rexp,abm,abp,omega,n,sefmnew,sefpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION rexp(n/2)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm0
      DOUBLE PRECISION omega(n),abm(n),abp(n),beta
      DOUBLE PRECISION sefmnew(n),sefpnew(n)
      EXTERNAL cspec,ferm0
      
*abm and abp are the heavy boson spectral functions        

      DO i=1,n/2
         exom=omega(i)
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=abm(m)+(abm(m+1)-abm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=abm(l)+(abm(l+1)-abm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(omega(j+1))*ferm0(-beta*omega(j+1))*wert2+
     &              cspec(omega(j))*ferm0(-beta*omega(j))*wert1)
c               write(86,*) omega(j),wert2
            ELSE
cc            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=abm(m+k)
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(freq2)*ferm0(-beta*freq2)*wertn+
     &                 cspec(freq1)*ferm0(-beta*freq1)*wert1)
                  freq1=freq2
                  wert1=wertn
               END DO
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (cspec(omega(j+1))*ferm0(-beta*omega(j+1))*wert2+
     &              cspec(freq1)*ferm0(-beta*freq1)*wert1)
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sefmnew(i)=aux
      END DO
*******************positive externa frequencies**********************  
      DO i=n/2+1,n
         exom=omega(i)
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=abp(l)+(abp(l+1)-abp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(omega(j+1))*ferm0(beta*omega(j+1))*wert2+
     &              cspec(omega(j))*ferm0(beta*omega(j))*wert1)
c     write(86,*) omega(j),wert2
            ELSE
c               IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
c     write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=abp(m+k)
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(freq2)*ferm0(beta*freq2)*wertn+
     &                 cspec(freq1)*ferm0(beta*freq1)*wert1)
                  freq1=freq2
                  wert1=wertn
               END DO
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (cspec(omega(j+1))*ferm0(beta*omega(j+1))*wert2+
     &              cspec(freq1)*ferm0(beta*freq1)*wert1)
            ENDIF
 30         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
         END DO
         sefpnew(i)=aux
      END DO
      DO i=1,n/2
         sefpnew(i)=rexp(i)*sefmnew(i)
         sefmnew(i+n/2)=rexp(n/2-i+1)*sefpnew(i+n/2)
      END DO
      
      RETURN
      END
************************************************************************
*****************************SIGMAB_HEAVY_U*****************************
************************************************************************
      SUBROUTINE SIGMAU(beta,rexp,afm,afp,omega,n,sebmnew,sebpnew)      
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION rexp(n/2)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm0
      DOUBLE PRECISION omega(n),afm(n),afp(n),beta
      DOUBLE PRECISION sebmnew(n),sebpnew(n)
      EXTERNAL cspec,ferm0      
        

      DO i=1,n/2
         exom=omega(i)
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(-omega(j+1))*ferm0(-beta*omega(j+1))*wert2+
     &              cspec(-omega(j))*ferm0(-beta*omega(j))*wert1)
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afm(m+k)
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(-freq2)*ferm0(-beta*freq2)*wertn+
     &                 cspec(-freq1)*ferm0(-beta*freq1)*wert1)
                  freq1=freq2
                  wert1=wertn
               END DO
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (cspec(-omega(j+1))*ferm0(-beta*omega(j+1))*wert2+
     &              cspec(-freq1)*ferm0(-beta*freq1)*wert1)
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sebmnew(i)=aux
      END DO
*******************positive externa frequencies**********************  
      DO i=n/2+1,n
         exom=omega(i)
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(-omega(j+1))*ferm0(beta*omega(j+1))*wert2+
     &              cspec(-omega(j))*ferm0(beta*omega(j))*wert1)
c     write(86,*) omega(j),wert2
            ELSE
c               IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
c     write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(-freq2)*ferm0(beta*freq2)*wertn+
     &                 cspec(-freq1)*ferm0(beta*freq1)*wert1)
                  freq1=freq2
                  wert1=wertn
               END DO
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (cspec(-omega(j+1))*ferm0(beta*omega(j+1))*wert2+
     &              cspec(-freq1)*ferm0(beta*freq1)*wert1)
            ENDIF
 30         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
         END DO
         sebpnew(i)=aux
      END DO
      DO i=1,n/2
         sebpnew(i)=rexp(i)*sebmnew(i)
         sebmnew(i+n/2)=rexp(n/2-i+1)*sebpnew(i+n/2)
      END DO
      
      RETURN
      END

************************************************************************
************************************************************************
***************Ankopplung an bosonischen Spinbad************************
************************************************************************
      SUBROUTINE SIGMAUbath(beta,rexp,afm,afp,omega,n,sefm,sefp)      
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION rexp(n/2)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex,D,xxx,www
      DOUBLE PRECISION aphi,bose0
      DOUBLE PRECISION omega(n),afm(n),afp(n),beta
      DOUBLE PRECISION sefm(n),sefp(n)
      EXTERNAL aphi,bose0      
        
      D=omega(n)
      DO i=1,n/2
         exom=omega(i)
c         xxx=aphi(0.0d0)*afm(i)
c         www=xxx*D*(-1.0d0)!to be added to integral
         www=0.0d0
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (bose0(-omega(j+1))*wert2*aphi(omega(j+1))+
     &              bose0(-omega(j))*wert1*aphi(omega(j)))
**********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afm(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (bose0(-freq2)*wertn*aphi(freq2)+
     &                 bose0(-freq1)*wert1*aphi(freq1))
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (bose0(-omega(j+1))*wert2*aphi(omega(j+1))+
     &              bose0(-freq1)*wert1*aphi(freq1))
**********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sefm(i)=-(aux+www)
      END DO
*******************positive externa frequencies**********************  
      DO i=n/2+1,n
         exom=omega(i)
c         xxx=aphi(0.0d0)*afp(i)
c         www=xxx*D*(-1.0d0)!to be added to integral
         www=0.0d0
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (bose0(omega(j+1))*wert2*aphi(omega(j+1))+
     &              bose0(omega(j))*wert1*aphi(omega(j)))
**********************************************************integration
c     write(86,*) omega(j),wert2
            ELSE
c               IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
c     write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (bose0(freq2)*wertn*aphi(freq2)+
     &                 bose0(freq1)*wert1*aphi(freq1))
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (bose0(omega(j+1))*wert2*aphi(omega(j+1))+
     &              bose0(freq1)*wert1*aphi(freq1))
**********************************************************integration
            ENDIF
 30         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
         END DO
         sefp(i)=aux+www
      END DO
      DO i=1,n/2
         sefp(i)=rexp(i)*sefm(i)
         sefm(i+n/2)=rexp(n/2-i+1)*sefp(i+n/2)
      END DO
      
      RETURN
      END
