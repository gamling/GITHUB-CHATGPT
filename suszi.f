      SUBROUTINE SUSCEPTIBILITY(beta,afp,afm,omega,n,suszi)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,mm,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION omega(n),afp(n),beta
      DOUBLE PRECISION suszi(n),wfreq1,wfreq2
      DOUBLE PRECISION afm(n),wert1m,wert2m,wmfreq1,wmfreq2,wertnm
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
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
             wert1m=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
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
               wert2m=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               wert2m=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               wert2m=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (afm(j+1)*wert2-wert2m*afp(j+1)
     &          +afm(j)*wert1-wert1m*afp(j))
***********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               wfreq1=afp(j)
               wmfreq1=afm(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  CALL LOCATE(omega,n,freq2,mm)
                  wfreq2=afp(mm)+(afp(mm+1)-afp(mm))/
     &              (omega(mm+1)-omega(mm))
     &           *(freq2-omega(mm)) !interpolate A_f
                  wmfreq2=afm(mm)+(afm(mm+1)-afm(mm))/
     &              (omega(mm+1)-omega(mm))
     &           *(freq2-omega(mm)) !interpolate A_f
                  wertn=afp(m+k)
                  wertnm=afm(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (wmfreq2*wertn-wfreq2*wertnm
     &             +wmfreq1*wert1-wfreq1*wert1m)
***********************************************************integration
                  freq1=freq2
                  wfreq1=wfreq2
                  wert1=wertn
                  wmfreq1=wmfreq2
                  wert1m=wertnm
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (afm(j+1)*wert2-wert2m*afp(j+1)
     &              +wmfreq1*wert1-wfreq1*wert1m)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
            wert1m=wert2m
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         suszi(i)=aux
      END DO
    
      RETURN
      END
