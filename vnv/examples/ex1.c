/******************************************************************************
 * Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
 * HYPRE Project Developers. See the top-level COPYRIGHT file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 ******************************************************************************/

/*
   Example 1

   Interface:    Structured interface (Struct)

   Compile with: make ex1 (may need to edit HYPRE_DIR in Makefile)

   Sample run:   mpirun -np 2 ex1

   Description:  This is a two processor example.  Each processor owns one
                 box in the grid.  For reference, the two grid boxes are those
                 in the example diagram in the struct interface chapter
                 of the User's Manual. Note that in this example code, we have
                 used the two boxes shown in the diagram as belonging
                 to processor 0 (and given one box to each processor). The
                 solver is PCG with no preconditioner.

                 We recommend viewing examples 1-4 sequentially for
                 a nice overview/tutorial of the struct interface.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Struct linear solvers header */
#include "HYPRE_struct_ls.h"
#include "ex.h"
#include "VnV.h"


/**
 * @title Hypre Example 1: The Struct Interface
 * @shortTitle Example 1
 *
 * This is a two processor example.  Each processor owns one
 * box in the grid.  For reference, the two grid boxes are those
 * in the example diagram in the struct interface chapter
 * of the User's Manual. Note that in this example code, we have
 * used the two boxes shown in the diagram as belonging
 * to processor 0 (and given one box to each processor). The
 * solver is PCG with no preconditioner.
 *

 *    
 */ 
INJECTION_EXECUTABLE(HYPRE_EX1)


/**
 * This comment gets parsed, but does not get used for anything 
 * yet. It could eventually be used to describe how the application
 * links with the subpackage I guess. 
 *
 * .. comment:: 
 *    
 *    This command tells the VnV Preprocessor that this executable 
 *    makes calls to a VnV package called VnVHypre. At runtime, VnV
 *    will call the registration function for this library and also
 *    for any subpacakges. Note, if you list a subpackage here, but 
 *    do not link a library that has that package, you will get an 
 *    error saying something along the lines of "undefined symbol 
 *    vnv_registration_<packagename>".   
 **/
 INJECTION_SUBPACKAGE(HYPRE_EX1, VnVHypre)

	
/**
 * @title Hypre Example 1
 * @shorttitle Hypre Example 1
 *
 * This is a two processor example.  Each processor owns one
 * box in the grid.  For reference, the two grid boxes are those
 * in the example diagram in the struct interface chapter
 * of the User's Manual. Note that in this example code, we have
 * used the two boxes shown in the diagram as belonging
 * to processor 0 (and given one box to each processor). The
 * solver is PCG with no preconditioner.
 * 
 * .. comment::
 *
 *    This comment shows up under the "packages" tab in the final 
 *    report. You can use this comment to print out any information
 *    you think is pertinent to the final report. For example, you
 *    could print out license information, version information, etc.
 *    C++ users have access to an engine* object that they can use to 
 *    write out variables to be used in this comment just like in the 
 *    vnv tests. 
 *
 */ 
 INJECTION_OPTIONS(HYPRE_EX1, "{ \"type\" : \"object\" }", void){
  
        // The options function gets passed a json-like object that 
	// is pre-validated against the schema passed above. In this 
	// case the schema is just a generic object, so all config is 
	// accepted. Hypre is compiled in C, so in this case the object is a "cjson"
	// object called json. The cjson api is described here (TODO -- document CJson.h)
	// TODO -Option to force C API in cases where we compile with C++ compiler. 

 	//This one just dumps to screen as a demo.	 

	// You can return a pointer to any object you like here. VnV
	// will store it in static space. You can fetch it at any time
	// using the INJECTION_GET_OPTIONS_OBJECT(HYPRE_EX1) command.
	
	return NULL;
}




int main (int argc, char *argv[])
{
   int  i, j, myid, num_procs;

   HYPRE_Int vis = 0;

   HYPRE_StructGrid     grid;
   HYPRE_StructStencil  stencil;
   HYPRE_StructMatrix   A;
   HYPRE_StructVector   b;
   HYPRE_StructVector   x;
   HYPRE_StructSolver   solver;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


   /* Initialize HYPRE */
   HYPRE_Init();

   /* Print GPU info */
   /* HYPRE_PrintDeviceInfo(); */

   HYPRE_Int vnv_index = -1;
   /* Parse command line */
   {
      HYPRE_Int arg_index = 0;
      HYPRE_Int print_usage = 0;

      while (arg_index < argc)
      {
         if ( strcmp(argv[arg_index], "-vis") == 0 )
         {
            arg_index++;
            vis = 1;
         }
         else if ( strcmp(argv[arg_index], "-help") == 0 )
         {
            print_usage = 1;
            break;
         }
	 else if (strcmp(argv[arg_index], "-vnv") == 0 )
	 {
	    vnv_index = ++arg_index;
	    ++arg_index;
	 }
	 else
         {
            arg_index++;
         }
      }

      if ((print_usage) && (myid == 0))
      {
         printf("\n");
         printf("Usage: %s [<options>]\n", argv[0]);
         printf("\n");
         printf("  -vis : save the solution for GLVis visualization\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }
   }

   if (vnv_index > -1 ) {
	INJECTION_INITIALIZE(HYPRE_EX1, &argc, &argv, argv[vnv_index]);   
   }

   if (num_procs != 2)
   {
      if (myid == 0) { printf("Must run with 2 processors!\n"); }
   } 
   else {
   /* 1. Set up a grid. Each processor describes the piece
      of the grid that it owns. */
   {
      /* Create an empty 2D grid object */
      HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

      /* Add boxes to the grid */
      if (myid == 0)
      {
         HYPRE_Int ilower[2] = {-3, 1}, iupper[2] = {-1, 2};
         HYPRE_StructGridSetExtents(grid, ilower, iupper);
      }
      else if (myid == 1)
      {
         HYPRE_Int ilower[2] = {0, 1}, iupper[2] = {2, 4};
         HYPRE_StructGridSetExtents(grid, ilower, iupper);
      }

      /* This is a collective call finalizing the grid assembly.
         The grid is now ``ready to be used'' */
      HYPRE_StructGridAssemble(grid);
   }

   /* 2. Define the discretization stencil */
   {
      /* Create an empty 2D, 5-pt stencil object */
      HYPRE_StructStencilCreate(2, 5, &stencil);

      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */
      {
         HYPRE_Int entry;
         HYPRE_Int offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};

         /* Assign each of the 5 stencil entries */
         for (entry = 0; entry < 5; entry++)
         {
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
         }
      }
   }

   /* 3. Set up a Struct Matrix */
   {
      /* Create an empty matrix object */
      HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);

      /* Indicate that the matrix coefficients are ready to be set */
      HYPRE_StructMatrixInitialize(A);

      /* Set the matrix coefficients.  Each processor assigns coefficients
         for the boxes in the grid that it owns. Note that the coefficients
         associated with each stencil entry may vary from grid point to grid
         point if desired.  Here, we first set the same stencil entries for
         each grid point.  Then we make modifications to grid points near
         the boundary. */
      if (myid == 0)
      {
         HYPRE_Int ilower[2] = {-3, 1}, iupper[2] = {-1, 2};
         HYPRE_Int stencil_indices[5] = {0, 1, 2, 3, 4}; /* labels for the stencil entries -
                                                  these correspond to the offsets
                                                  defined above */
         HYPRE_Int nentries = 5;
         HYPRE_Int nvalues  = 30; /* 6 grid points, each with 5 stencil entries */
         /* double values[30]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(30 * sizeof(double));

         /* We have 6 grid points, each with 5 stencil entries */
         for (i = 0; i < nvalues; i += nentries)
         {
            values[i] = 4.0;
            for (j = 1; j < nentries; j++)
            {
               values[i + j] = -1.0;
            }
         }

         HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                        stencil_indices, values);

         free(values);
      }
      else if (myid == 1)
      {
         HYPRE_Int ilower[2] = {0, 1}, iupper[2] = {2, 4};
         HYPRE_Int stencil_indices[5] = {0, 1, 2, 3, 4};
         HYPRE_Int nentries = 5;
         HYPRE_Int nvalues  = 60; /* 12 grid points, each with 5 stencil entries */
         /* double values[60]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(60 * sizeof(double));

         for (i = 0; i < nvalues; i += nentries)
         {
            values[i] = 4.0;
            for (j = 1; j < nentries; j++)
            {
               values[i + j] = -1.0;
            }
         }

         HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                        stencil_indices, values);

         free(values);
      }

      /* Set the coefficients reaching outside of the boundary to 0 */
      if (myid == 0)
      {
         /* double values[3]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(3 * sizeof(double));
         for (i = 0; i < 3; i++)
         {
            values[i] = 0.0;
         }
         {
            /* values below our box */
            HYPRE_Int ilower[2] = {-3, 1}, iupper[2] = {-1, 1};
            HYPRE_Int stencil_indices[1] = {3};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         {
            /* values to the left of our box */
            HYPRE_Int ilower[2] = {-3, 1}, iupper[2] = {-3, 2};
            HYPRE_Int stencil_indices[1] = {1};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         {
            /* values above our box */
            HYPRE_Int ilower[2] = {-3, 2}, iupper[2] = {-1, 2};
            HYPRE_Int stencil_indices[1] = {4};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         free(values);
      }
      else if (myid == 1)
      {
         /* double values[4]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(4 * sizeof(double));
         for (i = 0; i < 4; i++)
         {
            values[i] = 0.0;
         }
         {
            /* values below our box */
            HYPRE_Int ilower[2] = {0, 1}, iupper[2] = {2, 1};
            HYPRE_Int stencil_indices[1] = {3};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         {
            /* values to the right of our box */
            HYPRE_Int ilower[2] = {2, 1}, iupper[2] = {2, 4};
            HYPRE_Int stencil_indices[1] = {2};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         {
            /* values above our box */
            HYPRE_Int ilower[2] = {0, 4}, iupper[2] = {2, 4};
            HYPRE_Int stencil_indices[1] = {4};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         {
            /* values to the left of our box
               (that do not border the other box on proc. 0) */
            HYPRE_Int ilower[2] = {0, 3}, iupper[2] = {0, 4};
            HYPRE_Int stencil_indices[1] = {1};
            HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1,
                                           stencil_indices, values);
         }
         free(values);
      }

      /* This is a collective call finalizing the matrix assembly.
         The matrix is now ``ready to be used'' */
      HYPRE_StructMatrixAssemble(A);
   }

   /* 4. Set up Struct Vectors for b and x.  Each processor sets the vectors
      corresponding to its boxes. */
   {
      /* Create an empty vector object */
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);

      /* Set the vector coefficients */
      if (myid == 0)
      {
         HYPRE_Int ilower[2] = {-3, 1}, iupper[2] = {-1, 2};
         /* double values[6]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(6 * sizeof(double)); /* 6 grid points */

         for (i = 0; i < 6; i ++)
         {
            values[i] = 1.0;
         }
         HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

         for (i = 0; i < 6; i ++)
         {
            values[i] = 0.0;
         }
         HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
         free(values);
      }
      else if (myid == 1)
      {
         HYPRE_Int ilower[2] = {0, 1}, iupper[2] = {2, 4};
         /* double values[12]; OK to use constant-length arrays for CPUs */
         double *values = (double *) malloc(12 * sizeof(double)); /* 12 grid points */

         for (i = 0; i < 12; i ++)
         {
            values[i] = 1.0;
         }
         HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

         for (i = 0; i < 12; i ++)
         {
            values[i] = 0.0;
         }
         HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
         free(values);
      }

      /* This is a collective call finalizing the vector assembly.
         The vectors are now ``ready to be used'' */
      HYPRE_StructVectorAssemble(b);
      HYPRE_StructVectorAssemble(x);
   }

   /* 5. Set up and use a solver (See the Reference Manual for descriptions
      of all of the options.) */
   {
      /* Create an empty PCG Struct solver */
      HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);

      /* Set some parameters */
      HYPRE_StructPCGSetTol(solver, 1.0e-06); /* convergence tolerance */
      HYPRE_StructPCGSetPrintLevel(solver, 2); /* amount of info. printed */

      /* Setup and solve */
      HYPRE_StructPCGSetup(solver, A, b, x);
      HYPRE_StructPCGSolve(solver, A, b, x);
   }

   /* Save the solution for GLVis visualization, see vis/glvis-ex1.sh */
   if (vis)
   {
#ifdef HYPRE_EXVIS
      GLVis_PrintStructGrid(grid, "vis/ex1.mesh", myid, NULL, NULL);
      GLVis_PrintStructVector(x, "vis/ex1.sol", myid);
      GLVis_PrintData("vis/ex1.data", myid, num_procs);
#endif
   }

   /* Free memory */
   HYPRE_StructGridDestroy(grid);
   HYPRE_StructStencilDestroy(stencil);
   HYPRE_StructMatrixDestroy(A);
   HYPRE_StructVectorDestroy(b);
   HYPRE_StructVectorDestroy(x);
   HYPRE_StructPCGDestroy(solver);
  }

   //FINALIZE VNV
   if (vnv_index > -1 ) {
	INJECTION_FINALIZE(HYPRE_EX1);
   }

   /* Finalize Hypre */
   HYPRE_Finalize();
   

   /* Finalize MPI */
   MPI_Finalize();



   return (0);
}
