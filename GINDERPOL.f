C#################################################################
      SUBROUTINE GINDERPOL(grid0,n0,rneugrid,n2,func0,rneufunc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
       parameter (pi=3.1415926535898d0) 
       
       DIMENSION grid0(n0),func0(n0)
       DIMENSION grid(n0),rneugrid(n2),func(n0),rneufunc(n2)

C     make sure no points on grid coincide
       n1=1
       j=1
       grid(1)=grid0(1)
       func(1)=func0(1)
       DO i=2,n0
          j=j+1
          IF((grid0(i)-grid0(i-1)).LT.1.0d-9) then
          j=j-1
          else
             grid(j)=grid0(i)
             func(j)=func0(i)
             n1=n1+1
          endif
       END DO

       DO i=1,n2
          x=rneugrid(i)
          CALL LOCATE(grid,n1,x,j)
          
          IF(j.eq.0) THEN
             rneufunc(i)=0.0d0
          ELSEIF(j.eq.1) THEN
             rneufunc(i)=(func(2)-func(1))/(grid(2)-grid(1))*(x-grid(1))
          ELSEIF(j.eq.n1-1) THEN
             rneufunc(i)=(func(n1)-func(n1-1))/(grid(n1)-grid(n1-1))
     &            *(x-grid(n1-1))
           ELSEIF(j.eq.n1) THEN
              rneufunc(i)=0.0d0
           ELSE
              IF(grid(j+1).EQ.rneugrid(i)) THEN
                 rneufunc(i)=func(j+1)
              ELSE
                 d21=grid(j)-grid(j-1)
                 d32=grid(j+1)-grid(j)
                 d31=grid(j+1)-grid(j-1)
                 x1=x-grid(j-1)
                 x2=x-grid(j)
                 x3=x-grid(j+1)
                 rneufunc(i)=x2*x3/d21/d31*func(j-1)-
     &                x1*x3/d21/d32*func(j)+
     &                x1*x2/d31/d32*func(j+1)
              ENDIF
           ENDIF
        END DO

       RETURN
       END

      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
