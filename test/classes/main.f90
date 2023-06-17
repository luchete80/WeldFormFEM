  program main

    integer*8 myshape

    call readstep('cube.stp'//CHAR(0),myshape)

    call pass_it(myshape)

  end

  subroutine hello
    write (*,*) 'Hello from FORTRAN 77!'
  end subroutine