void atb_par(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
{
  int i, j, k;

  if( Ni == 1 && Nj == 1 && Nk == 1 ) {
    return B[0]*C[0];
  }
  else if( Ni*Nj <= 64*64 && Nk >= 2048 ) {

    // run k in parallel

  }
  else if( Ni <= 16 && Nj >= 8192 ) {

    // tile j in parallel

  }
  else if( Ni <= 8 && Nj >= 512 ) {

    // tile j in parallel

  }
  else {

    if( Ni <= 4 ) {

      if( Nk <= 2 ) {
        // no tiling no unrolling
      }
      else {
        // no tiling, k unrolled by 2

      }


    }
    else if( Nk < 256 ) {
      // i tiling with unrolling

      if( Nk <= 2 ) {

        // i tiling without unrolling

      }
      else {

        // i tiling with unrolling

      }
    }
    else {
      // ik tiling with unrolling

      if( Nk <= 2 ) {

        // ik tiling without unrolling

      }
      else {

        // ik tiling with unrolling

      }

    }


  }
}

