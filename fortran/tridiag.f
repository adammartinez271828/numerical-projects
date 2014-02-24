      SUBROUTINE tridiag(a,b,c,r,u,n)
      implicit none
      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet,gam(n)

c     a is AL, b is AM, c is AR
c     r is rhs, u is output, n is dimension

      if(b(1).eq.0.)pause 'tridiag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +%6V+j)D2.

