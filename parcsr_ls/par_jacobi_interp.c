/*BHEADER**********************************************************************
 * (c) 2006   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/

#include "headers.h"
#include "par_amg.h"

/* #define HYPRE_JACINT_PRINT_ROW_SUMS  */
/* #define HYPRE_JACINT_PRINT_SOME_ROWS */
/* #define HYPRE_JACINT_PRINT_MATRICES  */
#define HYPRE_JACINT_PRINT_DIAGS

void hypre_BoomerAMGJacobiInterp( hypre_ParCSRMatrix * A,
                                  hypre_ParCSRMatrix ** P,
                                  hypre_ParCSRMatrix * S,
                                  int * CF_marker, int level,
                                  double truncation_threshold,
                                  double truncation_threshold_minus )
/* nji steps of Jacobi interpolation, with nji presently just set in the code.*/
{
   double weight_AF = 1.0;  /* weight multiplied by A's fine row elements */
   double weight_AFF = 1.01;  /* weight multiplied by A's fine row + fine column elements */
   int nji = 1;
   int iji;

   for ( iji=0; iji<nji; ++iji )
   {
      hypre_BoomerAMGJacobiInterp_1( A, P, S, CF_marker, level,
                                     truncation_threshold, truncation_threshold_minus,
                                     weight_AF, weight_AFF );
   }

}

void hypre_BoomerAMGJacobiInterp_1( hypre_ParCSRMatrix * A,
                                    hypre_ParCSRMatrix ** P,
                                    hypre_ParCSRMatrix * S,
                                    int * CF_marker, int level,
                                    double truncation_threshold,
                                    double truncation_threshold_minus,
                                    double weight_AF, double weight_AFF
   )
/* One step of Jacobin interpolation.
   A is the linear system.
   P is an interpolation matrix, input and output
   CF_marker identifies coarse and fine points
   If we imagine P and A as split into coarse and fine submatrices,

       [ AFF  AFC ]   [ AF ]            [ IFC ]
   A = [          ] = [    ] ,      P = [     ]
       [ ACF  ACC ]   [ AC ]            [ ICC ]
   (note that ICC is an identity matrix, applied to coarse points only)
   then this function computes

   IFCnew = IFCold - DFF(-1) * ( AFF*IFCold + AFC )
          = IFCold - DFF(-1) * AF * Pold)
   where DFF is the diagonal of AFF, (-1) represents the inverse, and
   where "old" denotes a value on entry to this function, "new" a returned value.

 */
{
   hypre_ParCSRMatrix * Pnew;
   hypre_ParCSRMatrix * C;
   hypre_CSRMatrix *P_diag = hypre_ParCSRMatrixDiag(*P);
   hypre_CSRMatrix *P_offd = hypre_ParCSRMatrixOffd(*P);
   double          *P_diag_data = hypre_CSRMatrixData(P_diag);
   int             *P_diag_i = hypre_CSRMatrixI(P_diag);
   int             *P_diag_j = hypre_CSRMatrixJ(P_diag);
   double          *P_offd_data = hypre_CSRMatrixData(P_offd);
   int             *P_offd_i = hypre_CSRMatrixI(P_offd);
   hypre_CSRMatrix *C_diag;
   hypre_CSRMatrix *C_offd;
   hypre_CSRMatrix *Pnew_diag;
   hypre_CSRMatrix *Pnew_offd;
   int	num_rows_diag_P = hypre_CSRMatrixNumRows(P_diag);
   int i;
   int Jnochanges=0, Jchanges, Pnew_num_nonzeros;
   int CF_coarse=0;
   int * J_marker = hypre_CTAlloc( int, num_rows_diag_P );
   int nc, ncmax, ncmin, nc1;
   int num_procs, my_id;
   MPI_Comm comm = hypre_ParCSRMatrixComm( A );
#ifdef HYPRE_JACINT_PRINT_ROW_SUMS
   int m, nmav, npav;
   double PIi, PIimax, PIimin, PIimav, PIipav, randthresh;
   double eps = 1.0e-17;
#endif
#ifdef HYPRE_JACINT_PRINT_MATRICES
   char filename[80];
#endif
#ifdef HYPRE_JACINT_PRINT_SOME_ROWS
   int sample_rows[50], n_sample_rows=0, isamp;
#endif

   MPI_Comm_size(comm, &num_procs);   
   MPI_Comm_rank(comm,&my_id);


   for ( i=0; i<num_rows_diag_P; ++i )
   {
      J_marker[i] = CF_marker[i];
      if (CF_marker[i]>=0) ++CF_coarse;
   }
#ifdef HYPRE_JACINT_PRINT_DIAGS
   printf("%i %i Jacobi_Interp_1, P has %i+%i=%i nonzeros\n", my_id, level,
          hypre_CSRMatrixNumNonzeros(P_diag), hypre_CSRMatrixNumNonzeros(P_offd),
          hypre_CSRMatrixNumNonzeros(P_diag)+hypre_CSRMatrixNumNonzeros(P_offd)  );
#endif

   /* row sum computations, for output */
#ifdef HYPRE_JACINT_PRINT_ROW_SUMS
   PIimax=-1.0e12, PIimin=1.0e12, PIimav=0, PIipav=0;
   nmav=0, npav=0;
   for ( i=0; i<num_rows_diag_P; ++i )
   {
      PIi = 0;  /* i-th value of P*1, i.e. sum of row i of P */
      for ( m=P_diag_i[i]; m<P_diag_i[i+1]; ++m )
         PIi += P_diag_data[m];
      for ( m=P_offd_i[i]; m<P_offd_i[i+1]; ++m )
         PIi += P_offd_data[m];
      if (CF_marker[i]<0)
      {
         PIimax = hypre_max( PIimax, PIi );
         PIimin = hypre_min( PIimin, PIi );
         if (PIi<=1-eps) { PIimav+=PIi; ++nmav; };
         if (PIi>=1+eps) { PIipav+=PIi; ++npav; };
      }
   }
   if ( nmav>0 ) PIimav = PIimav/nmav;
   if ( npav>0 ) PIipav = PIipav/npav;
   printf("%i %i P in max,min row sums %e %e\n", my_id, level, PIimax, PIimin );
#endif

   ncmax=0; ncmin=num_rows_diag_P; nc1=0;
   for ( i=0; i<num_rows_diag_P; ++i )
      if (CF_marker[i]<0)
      {
         nc = P_diag_i[i+1] - P_diag_i[i];
         if (nc<=1)
         {
            ++nc1;
         }
         ncmax = hypre_max( nc, ncmax );
         ncmin = hypre_min( nc, ncmin );
      }
#if 0
   /* a very agressive reduction in how much the Jacobi step does: */
   for ( i=0; i<num_rows_diag_P; ++i )
      if (CF_marker[i]<0)
      {
         nc = P_diag_i[i+1] - P_diag_i[i];
         if (nc>ncmin+1)
            /*if ( nc > ncmin + 0.5*(ncmax-ncmin) )*/
         {
            J_marker[i] = 1;
            ++Jnochanges;
         }
      }
#endif
   Jchanges = num_rows_diag_P - Jnochanges - CF_coarse;

#ifdef HYPRE_JACINT_PRINT_SOME_ROWS
   printf("some rows to be changed: ");
   randthresh = 15/(double)Jchanges;
   for ( i=0; i<num_rows_diag_P; ++i )
   {
      if ( J_marker[i]<0 )
      {
         if ( ((double)rand())/RAND_MAX < randthresh )
         {
            printf( "%i: ", i );
            for ( m=P_diag_i[i]; m<P_diag_i[i+1]; ++m )
               printf( " %i %f, ", P_diag_j[m], P_diag_data[m] );
            printf(";  ");
            sample_rows[n_sample_rows] = i;
            ++n_sample_rows;
         }
      }
   }
   printf("\n");
#endif
#ifdef HYPRE_JACINT_PRINT_DIAGS
   printf("%i %i P has %i rows, %i changeable, %i don't change-good, %i coarse\n",
          my_id, level, num_rows_diag_P, Jchanges, Jnochanges, CF_coarse );
   printf("%i %i min,max cols per row: %i, %i;  no.rows w.<=1 col: %i\n", my_id, level, ncmin, ncmax, nc1 );
#endif
#ifdef HYPRE_JACINT_PRINT_MATRICES
   if ( num_rows_diag_P <= 100 )
   {
      sprintf( filename, "Ain%i", level );
      hypre_ParCSRMatrixPrintIJ( A,0,0,filename);
      sprintf( filename, "Sin%i", level );
      hypre_ParCSRMatrixPrintIJ( S,0,0,filename);
      sprintf( filename, "Pin%i", level );
      hypre_ParCSRMatrixPrintIJ( *P,0,0,filename);
   }
#endif

   C = hypre_ParMatmul_FC( A, *P, J_marker, weight_AFF );
   /* hypre_parMatmul_FC creates and returns C, a variation of the
      matrix product A*P in which only the "Fine"-designated rows have
      been computed.  (all columns are Coarse because all columns of P
      are).  "Fine" is defined solely by the marker array, and for
      example could be a proper subset of the fine points of a
      multigrid hierarchy.
      As a matrix, C is the size of A*P.  But only the marked rows have
      been computed.
   */
#ifdef HYPRE_JACINT_PRINT_MATRICES
   sprintf( filename, "C%i", level );
   if ( num_rows_diag_P <= 100 ) hypre_ParCSRMatrixPrintIJ( C,0,0,filename);
#endif
   C_diag = hypre_ParCSRMatrixDiag(C);
   C_offd = hypre_ParCSRMatrixOffd(C);
#ifdef HYPRE_JACINT_PRINT_DIAGS
   printf("%i %i Jacobi_Interp_1 after matmul, C has %i+%i=%i nonzeros\n",
          my_id, level, hypre_CSRMatrixNumNonzeros(C_diag),
          hypre_CSRMatrixNumNonzeros(C_offd),
          hypre_CSRMatrixNumNonzeros(C_diag)+hypre_CSRMatrixNumNonzeros(C_offd)  );
#endif

   hypre_ParMatScaleDiagInv_F( C, A, weight_AF, J_marker );
   /* hypre_ParMatScaleDiagInv scales of its first argument by premultiplying with
      a submatrix of the inverse of the diagonal of its second argument.
      The marker array determines which diagonal elements are used.  The marker
      array should select exactly the right number of diagonal elements (the number
      of rows of AP_FC).
   */
#ifdef HYPRE_JACINT_PRINT_MATRICES
   sprintf( filename, "Cout%i", level );
   if ( num_rows_diag_P <= 100 )  hypre_ParCSRMatrixPrintIJ( C,0,0,filename);
#endif

   Pnew = hypre_ParMatMinus_F( *P, C, J_marker );
   /* hypre_ParMatMinus_F subtracts rows of its second argument from selected rows
      of its first argument.  The marker array determines which rows of the first
      argument are affected, and they should exactly correspond to all the rows
      of the second argument.
   */
   Pnew_diag = hypre_ParCSRMatrixDiag(Pnew);
   Pnew_offd = hypre_ParCSRMatrixOffd(Pnew);
   Pnew_num_nonzeros = hypre_CSRMatrixNumNonzeros(Pnew_diag)+hypre_CSRMatrixNumNonzeros(Pnew_offd);
#ifdef HYPRE_JACINT_PRINT_DIAGS
   printf("%i %i Jacobi_Interp_1 after MatMinus, Pnew has %i+%i=%i nonzeros\n",
          my_id, level, hypre_CSRMatrixNumNonzeros(Pnew_diag),
          hypre_CSRMatrixNumNonzeros(Pnew_offd), Pnew_num_nonzeros );
#endif

   hypre_ParCSRMatrixDestroy( C );
   /*>>> Memory leak - we can't destroy the old P, it still contains something
     >>> being used elsewhere, probably row_starts or something like it which
     >>> isn't "owned" by the new arrays...
     hypre_ParCSRMatrixDestroy( *P );*/

   /* Note that I'm truncating all the fine rows, not just the J-marked ones. */
   if ( Pnew_num_nonzeros < 10000 )
   {  /* ad-hoc attempt to reduce zero-matrix problems seen in testing..*/
      truncation_threshold = 1.0e-6 * truncation_threshold; 
      truncation_threshold_minus = 1.0e-6 * truncation_threshold_minus;
  }
   hypre_BoomerAMGTruncateInterp( Pnew, truncation_threshold,
                                  truncation_threshold_minus, CF_marker );
   *P = Pnew;

   P_diag = hypre_ParCSRMatrixDiag(*P);
   P_offd = hypre_ParCSRMatrixOffd(*P);
   P_diag_data = hypre_CSRMatrixData(P_diag);
   P_diag_i = hypre_CSRMatrixI(P_diag);
   P_diag_j = hypre_CSRMatrixJ(P_diag);
   P_offd_data = hypre_CSRMatrixData(P_offd);
   P_offd_i = hypre_CSRMatrixI(P_offd);

   /* row sum computations, for output */
#ifdef HYPRE_JACINT_PRINT_ROW_SUMS
   PIimax=-1.0e12, PIimin=1.0e12, PIimav=0, PIipav=0;
   nmav=0, npav=0;
   for ( i=0; i<num_rows_diag_P; ++i )
   {
      PIi = 0;  /* i-th value of P*1, i.e. sum of row i of P */
      for ( m=P_diag_i[i]; m<P_diag_i[i+1]; ++m )
         PIi += P_diag_data[m];
      for ( m=P_offd_i[i]; m<P_offd_i[i+1]; ++m )
         PIi += P_offd_data[m];
      if (CF_marker[i]<0)
      {
         PIimax = hypre_max( PIimax, PIi );
         PIimin = hypre_min( PIimin, PIi );
         if (PIi<=1-eps) { PIimav+=PIi; ++nmav; };
         if (PIi>=1+eps) { PIipav+=PIi; ++npav; };
      }
   }
   if ( nmav>0 ) PIimav = PIimav/nmav;
   if ( npav>0 ) PIipav = PIipav/npav;
   printf("%i %i P out max,min row sums %e %e\n", my_id, level, PIimax, PIimin );
#endif

#ifdef HYPRE_JACINT_PRINT_SOME_ROWS
   printf("some changed rows: ");
   for ( isamp=0; isamp<n_sample_rows; ++isamp )
   {
      i = sample_rows[isamp];
      printf( "%i: ", i );
      for ( m=P_diag_i[i]; m<P_diag_i[i+1]; ++m )
         printf( " %i %f, ", P_diag_j[m], P_diag_data[m] );
      printf(";  ");
   }
   printf("\n");
#endif
   ncmax=0; ncmin=num_rows_diag_P; nc1=0;
   for ( i=0; i<num_rows_diag_P; ++i )
      if (CF_marker[i]<0)
      {
         nc = P_diag_i[i+1] - P_diag_i[i];
         if (nc<=1) ++nc1;
         ncmax = hypre_max( nc, ncmax );
         ncmin = hypre_min( nc, ncmin );
      }
#ifdef HYPRE_JACINT_PRINT_DIAGS
   printf("%i %i P has %i rows, %i changeable, %i too good, %i coarse\n",
          my_id, level, num_rows_diag_P, num_rows_diag_P-Jnochanges-CF_coarse, Jnochanges, CF_coarse );
   printf("%i %i min,max cols per row: %i, %i;  no.rows w.<=1 col: %i\n", my_id, level, ncmin, ncmax, nc1 );

   printf("%i %i Jacobi_Interp_1 after truncation (%e), Pnew has %i+%i=%i nonzeros\n",
          my_id, level, truncation_threshold,
          hypre_CSRMatrixNumNonzeros(Pnew_diag), hypre_CSRMatrixNumNonzeros(Pnew_offd),
          hypre_CSRMatrixNumNonzeros(Pnew_diag)+hypre_CSRMatrixNumNonzeros(Pnew_offd)  );
#endif

   /* Programming Notes:
      1. Judging by around line 299 of par_interp.c, they typical use of CF_marker
      is that CF_marker>=0 means Coarse, CF_marker<0 means Fine.
   */
#ifdef HYPRE_JACINT_PRINT_MATRICES
   sprintf( filename, "Pout%i", level );
   if ( num_rows_diag_P <= 100 )  hypre_ParCSRMatrixPrintIJ( *P,0,0,filename);
#endif
      
}

void hypre_BoomerAMGTruncateInterp( hypre_ParCSRMatrix *P,
                                    double eps, double dlt,
                                    int * CF_marker )
/* Truncate the interpolation matrix P, but only in rows for which the
   marker is <0.  Truncation means that an element P(i,j) is set to 0 if
   P(i,j)>0 and P(i,j)<eps*max( P(i,j) )  or if
   P(i,j)>0 and P(i,j)<dlt*max( -P(i,j) )  or if
   P(i,j)<0 and P(i,j)>dlt*min( -P(i,j) )  or if
   P(i,j)<0 and P(i,j)>eps*min( P(i,j) )
      ( 0<eps,dlt<1, typically 0.1=dlt<eps=0.2, )
   The min and max are only computed locally, as I'm guessing that there isn't
   usually much to be gained (in the way of improved performance) by getting
   them perfectly right.
*/

/* The function hypre_BoomerAMGInterpTruncation in par_interp.c is
   very similar.  It looks at fabs(value) rather than separately
   dealing with value<0 and value>0 as recommended by Klaus Stuben,
   thus as this function does.  In this function, only "marked" rows
   are affected.  Lastly, in hypre_BoomerAMGInterpTruncation, if any
   element gets discarded, it reallocates arrays to the new size.
*/
{
   hypre_CSRMatrix *P_diag = hypre_ParCSRMatrixDiag(P);
   hypre_CSRMatrix *P_offd = hypre_ParCSRMatrixOffd(P);
   double          *P_diag_data = hypre_CSRMatrixData(P_diag);
   int             *P_diag_i = hypre_CSRMatrixI(P_diag);
   int             *P_diag_j = hypre_CSRMatrixJ(P_diag);
   double          *P_offd_data = hypre_CSRMatrixData(P_offd);
   int             *P_offd_i = hypre_CSRMatrixI(P_offd);
   int             *P_offd_j = hypre_CSRMatrixJ(P_offd);
   int             *new_P_diag_i;
   int             *new_P_offd_i;
   int	num_rows_diag_P = hypre_CSRMatrixNumRows(P_diag);
   int	num_rows_offd_P = hypre_CSRMatrixNumRows(P_offd);
   int num_nonzeros_diag = hypre_CSRMatrixNumNonzeros(P_diag);
   int num_nonzeros_offd = hypre_CSRMatrixNumNonzeros(P_offd);
   double vmax = 0.0;
   double vmin = 0.0;
   double v, old_sum, new_sum, scale, wmax, wmin;
   int i1, m, m1d, m1o;

   /* compute vmax = eps*max(P(i,j)), vmin = eps*min(P(i,j)) */
   for ( i1 = 0; i1 < num_rows_diag_P; i1++ )
   {
      for ( m=P_diag_i[i1]; m<P_diag_i[i1+1]; ++m )
         {
            v = P_diag_data[m];
            vmax = hypre_max( v, vmax );
            vmin = hypre_min( v, vmin );
         }
      for ( m=P_offd_i[i1]; m<P_offd_i[i1+1]; ++m )
         {
            v = P_offd_data[m];
            vmax = hypre_max( v, vmax );
            vmin = hypre_min( v, vmin );
         }
   }
   if ( vmax <= 0.0 ) vmax =  1.0;  /* make sure no v is v>vmax if no v is v>0 */
   if ( vmin >= 0.0 ) vmin = -1.0;  /* make sure no v is v<vmin if no v is v<0 */
   wmax = - dlt * vmin;
   wmin = - dlt * vmax;
   vmax *= eps;
   vmin *= eps;

   /* Repack the i,j,and data arrays so as to discard the small elements of P.
      Elements of Coarse rows (CF_marker>=0) are always kept.
      The arrays are not re-allocated, so there will generally be unused space
      at the ends of the arrays. */
   new_P_diag_i = hypre_CTAlloc( int, num_rows_diag_P+1 );
   new_P_offd_i = hypre_CTAlloc( int, num_rows_offd_P+1 );
   m1d = P_diag_i[0];
   m1o = P_offd_i[0];
   for ( i1 = 0; i1 < num_rows_diag_P; i1++ )
   {
      old_sum = 0;
      new_sum = 0;
      for ( m=P_diag_i[i1]; m<P_diag_i[i1+1]; ++m )
      {
         v = P_diag_data[m];
         old_sum += v;
         if ( CF_marker[i1]>=0 || ( v>=vmax && v>=wmax ) || ( v<=vmin && v<=wmin ) )
         {  /* keep v */
            new_sum += v;
            P_diag_j[m1d] = P_diag_j[m];
            P_diag_data[m1d] = P_diag_data[m];
            ++m1d;
         }
         else
         {  /* discard v */
            --num_nonzeros_diag;
         }
      }
      for ( m=P_offd_i[i1]; m<P_offd_i[i1+1]; ++m )
      {
         v = P_offd_data[m];
         old_sum += v;
         if ( CF_marker[i1]>=0 || ( v>=vmax && v>=wmax ) || ( v<=vmin && v<=wmin ) )
         {  /* keep v */
            new_sum += v;
            P_offd_j[m1o] = P_offd_j[m];
            P_offd_data[m1o] = P_offd_data[m];
            ++m1o;
         }
         else
         {  /* discard v */
            --num_nonzeros_offd;
         }
      }

      new_P_diag_i[i1+1] = m1d;
      if ( i1<num_rows_offd_P ) new_P_offd_i[i1+1] = m1o;

      /* rescale to keep row sum the same */
      if (new_sum!=0) scale = old_sum/new_sum; else scale = 1.0;
      for ( m=new_P_diag_i[i1]; m<new_P_diag_i[i1+1]; ++m )
         P_diag_data[m] *= scale;
      if ( i1<num_rows_offd_P ) /* this test fails when there is no offd block */
         for ( m=new_P_offd_i[i1]; m<new_P_offd_i[i1+1]; ++m )
            P_offd_data[m] *= scale;

   }

   for ( i1 = 1; i1 <= num_rows_diag_P; i1++ )
   {
      P_diag_i[i1] = new_P_diag_i[i1];
      if ( i1<=num_rows_offd_P && num_nonzeros_offd>0 ) P_offd_i[i1] = new_P_offd_i[i1];
   }
   hypre_TFree( new_P_diag_i );
   if ( num_rows_offd_P>0 ) hypre_TFree( new_P_offd_i );

   hypre_CSRMatrixNumNonzeros(P_diag) = num_nonzeros_diag;
   hypre_CSRMatrixNumNonzeros(P_offd) = num_nonzeros_offd;
   hypre_ParCSRMatrixSetDNumNonzeros( P );
   hypre_ParCSRMatrixSetNumNonzeros( P );

}