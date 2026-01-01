***********************************************************************
* Auxiliary Particle Bubble with defect propagators A^{\pm}           *
*Input omega,n,af(eps+omega),ab(eps)                                  *
*Output faltung                                                       *
***********************************************************************
      SUBROUTINE BUBBLE(omega,n,af,ab,faltung)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,mm,iscr(n)
      DOUBLE PRECISION omega(n),af(n),ab(n)
      DOUBLE PRECISION faltung(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION wfreq1,wfreq2,pi

      pi=4.0d0*DATAN(1.0d0)

      DO i=1,n
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
            wert1=af(m)+(af(m+1)-af(m))/(omega(m+1)-omega(m))
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
               wert2=af(l)+(af(l+1)-af(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (ab(j+1)*wert2+ab(j)*wert1)
***********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               wfreq1=ab(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  CALL LOCATE(omega,n,freq2,mm)
                  wfreq2=ab(mm)+(ab(mm+1)-ab(mm))/
     &              (omega(mm+1)-omega(mm))
     &           *(freq2-omega(mm)) !interpolate A_f
                  wertn=af(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (wfreq2*wertn+wfreq1*wert1)
***********************************************************integration
                  freq1=freq2
                  wfreq1=wfreq2
                  wert1=wertn
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (ab(j+1)*wert2+wfreq1*wert1)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         faltung(i)=aux/pi
      END DO
      
      RETURN
      END
