      program hw1
      implicit none
      double precision x,y,z
      
      write(6,*) 'Test 1'
      x = 1
      write(6,*) x
      x = 1.
      write(6,*) x
      x = 1.d0
      write(6,*) x
      write(6,*)
      
      write(6,*) 'Test 2'
      x = -3.d0
      write(6,*) x**2
      write(6,*) x**2.
      write(6,*) x**2.d0
      write(6,*) 
      
      write(6,*) 'Test 3'
      write(6,*) x**2.001d0
      write(6,*) x**2.001
      
      write(6,*) 'Test 4'
      write(6,*) x/0.d0
      
      stop
      end