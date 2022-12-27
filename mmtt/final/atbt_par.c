#include "lib_mylib_atbt.h"
#include "limits.h"

void atbt_par(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
{
    int i, j, k;

    // overflow detection
    int val0 = 999999;
    int val1 = 999999;

    int result0 = 0;
    int result1 = 0;

    if( Ni <= INT_MAX/Nj ) {

        val0 = Ni*Nj;
        result0 = 1;

        if( Nk <= INT_MAX/val0 ) {
            result1 = 1;
            val1 = val0*Nk;
        }

    }

    if( result1 && val1 <= 16*16*16 ) {

        // use sequential
        atbt_seq_ref( A, B, C, Ni, Nj, Nk );

    }
    else if( result0 && val0 <= 16*16 ) {

        if( val0 <= 8*8 ) {
             // use histogram code
             atbt_par_histogram( A, B, C, Ni, Nj, Nk );
        }
        else {

            if( Nk < 60 ) {

                // use serial code
                atbt_seq_ref( A, B, C, Ni, Nj, Nk );

            }
            else {

                // use k unrolled by 8
                atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

            }
            
        }

    }
    else if( Ni <= 64 ) {

        if( Nj <= 256 ) {

            // use k tiling
            if( Nk <= 2 ) {
                
                atbt_par_ref( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 4 ) {
                atbt_par_kunroll2( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 8 ) {
                atbt_par_kunroll4( A, B, C, Ni, Nj, Nk );

            }
            else {
                atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

            }

        }
        else if( Nj <= 512 ) {

            if( Nk <= 2048 ) {
                atbt_par_jtile16_kunroll8_ikj( A, B, C, Ni, Nj, Nk );

            }
            else {
                atbt_par_ktile_kunroll8( A, B, C, Ni, Nj, Nk );

            }

        }
        else if( Nk <= 1024 ) {

            if( Nk <= 2 ) {
                atbt_par_ref( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 4 ) {
                atbt_par_kunroll2( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 8 ) {
                atbt_par_kunroll4( A, B, C, Ni, Nj, Nk );

            }
            else {
                atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

            }

        }
        else {
            atbt_par_kjtile8_kunroll4( A, B, C, Ni, Nj, Nk );

        }

    }
    else if( Ni <= 128 ) {

        if( Nj <= 256 ) {

            // use k tiling
            if( Nk <= 2 ) {
                atbt_par_ref( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 4 ) {
                atbt_par_kunroll2( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 8 ) {
                atbt_par_kunroll4( A, B, C, Ni, Nj, Nk );

            }
            else {
                atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

            }

        }
        else if( Nj <= 512 ) {

            if( Nk <= 2 ) {
                atbt_par_ref( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 4 ) {
                atbt_par_kunroll2( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 8 ) {
                atbt_par_kunroll4( A, B, C, Ni, Nj, Nk );

            }
            else if( Nk <= 1024 ) {
                atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

            }
            else {
                atbt_par_ktile_kunroll8( A, B, C, Ni, Nj, Nk );

            }

        }

    }
    else if( Ni <= 256 ) {

        if( Nk <= 2 ) {
                atbt_par_ref( A, B, C, Ni, Nj, Nk );

        }
        else if( Nk <= 16 ) {
            atbt_par_kunroll4( A, B, C, Ni, Nj, Nk );

        }
        else {
            atbt_par_ktile_kunroll8( A, B, C, Ni, Nj, Nk );

        }

    }
    else if( Ni <= 1024 ) {

        if( Nj <= 64 ) {
            atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

        }
        else {
            atbt_par_ijtile16_kunroll4_kij( A, B, C, Ni, Nj, Nk );

        }

    }
    else if( Ni > 1024 && Nj > 1024 && Nk > 1024 ) {
        atbt_par_ikjtile16_kunroll4( A, B, C, Ni, Nj, Nk );

    }
    else if( Nj < 64 && Nk < 64 ) {
        atbt_par_kunroll8( A, B, C, Ni, Nj, Nk );

    }
    else {
        atbt_par_ijtile16_kunroll4_kij( A, B, C, Ni, Nj, Nk );

    }

}
