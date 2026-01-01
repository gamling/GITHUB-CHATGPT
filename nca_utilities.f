************************************************************************
************************************************************************
      subroutine ownsort(nges,n1,w,indx)
      implicit integer (i,j,k,l,m,n)
      implicit real*8 (a-h,o-z)
      dimension w(nges),aux(nges),indx(nges)
      
      
      n2=nges-n1
      l=1
      k=n1+1
      do i=1,nges
         mem=i
         if(w(l).lt.w(k)) then
            indx(i)=l
            aux(i)=w(l)
            if(l.eq.n1) then
               goto 30
            end if
            l=l+1
         else
            indx(i)=k
            aux(i)=w(k)
            if(k.eq.nges) then
               goto 40
            end if
            k=k+1
         end if
      end do

 30   continue
      do i=k,nges
         aux(i)=w(i)
         indx(i)=i
      end do
      goto 50
 40   continue
      do i=l,n1
         aux(i+n2)=w(i)
         indx(i+n2)=i
      end do

 50   continue
      do i=1,nges
         w(i)=aux(i)
      end do
      return
      end
************************************************************************
************************************************************************
*********************CALCULATES PARTICLE NUMBERS************************
      SUBROUTINE PN(exom,beta,n,omega,field,nf)
      IMPLICIT NONE
      INTEGER i,j,n,m,iscr(n),kk,k,l
      DOUBLE PRECISION exom,nf,beta,field(n),omega(n)
      DOUBLE PRECISION wert1,wert2,wertn,freq1,freq2
      DOUBLE PRECISION omex,aux
      DOUBLE PRECISION ferm0,pi
      PARAMETER(pi=3.1415926535898d0)
      EXTERNAL ferm0
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
            wert1=field(m)+(field(m+1)-field(m))/(omega(m+1)-omega(m))
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
               wert2=field(l)+(field(l+1)-field(l))/(omega(l+1)-omega(l)
     &              )*(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (ferm0(beta*omega(j+1))*wert2+
     &              ferm0(beta*omega(j))*wert1)
c     write(86,*) omega(j),wert2
            ELSE
c               IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
c     write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=field(m+k)
                  aux=aux+0.5d0*(freq2-freq1)*
     &                 (ferm0(beta*freq2)*wertn+
     &                 ferm0(beta*freq1)*wert1)
                  freq1=freq2
                  wert1=wertn
               END DO
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (ferm0(beta*omega(j+1))*wert2+
     &              ferm0(beta*freq1)*wert1)
            ENDIF
 30         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCpositive freq
         END DO
         nf=aux/pi
      RETURN
      END
