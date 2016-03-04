c compile it using: f95 -Wall -O3 test.f expm.o -lm -llapack -lblas -o test_expm
      program hello
         COMPLEX*16, allocatable::M(:,:)
         INTEGER DIM

         DIM = 3

         allocate(M(1:DIM,1:DIM))

         M(1,1) = 21
         M(1,2) = 17
         M(1,3) = 6

         M(2,1) = -5
         M(2,2) = -1
         M(2,3) = -6

         M(3,1) = 4
         M(3,2) = 4
         M(3,3) = 16

         call expm(M,DIM)
      end program hello
