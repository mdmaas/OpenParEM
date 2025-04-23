// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// This is a modified version of L2ZZErrorEstimator to 
// change the preconditioner, check for convergence,
// and return the error condition instead of the global error.

#include "OPEM_L2ZZErrorEstimator.hpp"

extern "C" PetscErrorCode hypre_ParCSRMatrixToMat(hypre_ParCSRMatrix *, Mat *, PetscInt, int, int, int);
extern "C" PetscErrorCode hypre_getSparseWidth (hypre_ParCSRMatrix *, PetscInt *);

bool OPEM_L2ZZErrorEstimator(BilinearFormIntegrator &flux_integrator,
                             const ParGridFunction &x,
                             ParFiniteElementSpace &smooth_flux_fes,
                             ParFiniteElementSpace &flux_fes,
                             Vector &errors_norm, 
                             double norm_p, double solver_tol, int solver_max_it,
                             double &finalResidualNorm, int debug_refine_preconditioner, int show_iteration_count)
{
   bool fail=false;

   // Compute fluxes in discontinuous space
   GridFunction flux(&flux_fes);
   flux = 0.0;

   ParFiniteElementSpace *xfes = x.ParFESpace();
   Array<int> xdofs, fdofs;
   Vector el_x, el_f;

   for (int i = 0; i < xfes->GetNE(); i++)
   {
      const DofTransformation* const xtrans = xfes->GetElementVDofs(i, xdofs);
      x.GetSubVector(xdofs, el_x);
      if (xtrans)
      {
         xtrans->InvTransformPrimal(el_x);
      }

      ElementTransformation *Transf = xfes->GetElementTransformation(i);
      flux_integrator.ComputeElementFlux(*xfes->GetFE(i), *Transf, el_x,
                                         *flux_fes.GetFE(i), el_f, false);

      const DofTransformation* const ftrans = flux_fes.GetElementVDofs(i, fdofs);
      if (ftrans)
      {
         ftrans->TransformPrimal(el_f);
      }
      flux.SetSubVector(fdofs, el_f);
   }

   // Assemble the linear system for L2 projection into the "smooth" space
   ParBilinearForm *a = new ParBilinearForm(&smooth_flux_fes);
   ParLinearForm *b = new ParLinearForm(&smooth_flux_fes);
   VectorGridFunctionCoefficient f(&flux);

   if (xfes->GetNE())
   {
      MFEM_VERIFY(smooth_flux_fes.GetFE(0) != NULL,
                  "Could not obtain FE of smooth flux space.");

      if (smooth_flux_fes.GetFE(0)->GetRangeType() == FiniteElement::SCALAR)
      {
         VectorMassIntegrator *vmass = new VectorMassIntegrator;
         vmass->SetVDim(smooth_flux_fes.GetVDim());
         a->AddDomainIntegrator(vmass);
         b->AddDomainIntegrator(new VectorDomainLFIntegrator(f));
      }
      else
      {
         a->AddDomainIntegrator(new VectorFEMassIntegrator);
         b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
      }
   }

   b->Assemble();
   a->Assemble();
   a->Finalize();

   // The destination of the projected discontinuous flux
   ParGridFunction smooth_flux(&smooth_flux_fes);
   smooth_flux = 0.0;

   HypreParMatrix* A = a->ParallelAssemble();
   HypreParVector* B = b->ParallelAssemble();
   HypreParVector* X = smooth_flux.ParallelProject();

   delete a;
   delete b;

   // solver
   HyprePCG *pcg = new HyprePCG(*A);
   pcg->SetTol(solver_tol);
   pcg->SetMaxIter(solver_max_it);
   pcg->SetPrintLevel(0);

   // preconditioner
   HypreDiagScale *diag=nullptr;
   HypreILU *ilu=nullptr;
   HypreBoomerAMG *boomer=nullptr;
   if (debug_refine_preconditioner == 0) {diag=new HypreDiagScale(*A); pcg->SetPreconditioner(*diag);}
   if (debug_refine_preconditioner == 1) {
      boomer=new HypreBoomerAMG(*A);
      boomer->SetRelaxType(6);
      boomer->SetPrintLevel(0);
      pcg->SetPreconditioner(*boomer);
   }

   // solve
   pcg->Mult(*B, *X);

   // check for convergence
   int numIterations;
   pcg->GetNumIterations(numIterations);
   if (numIterations >= solver_max_it) fail=true;

   pcg->GetFinalResidualNorm(finalResidualNorm);

   if (show_iteration_count == 1) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"            real part number of iterations: %d\n",numIterations);
   }

   if (show_iteration_count == 2) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"            imaginary part number of iterations: %d\n",numIterations);
   }

   // Extract the parallel grid function corresponding to the finite element
   // approximation X. This is the local solution on each processor.
   smooth_flux = *X;

   delete A;
   delete B;
   delete X;
   if (diag) {delete diag; diag=nullptr;}
   if (ilu) {delete ilu; ilu=nullptr;}
   if (boomer) {delete boomer; boomer=nullptr;}
   delete pcg;

   // Proceed through the elements one by one, and find the Lp norm differences
   // between the flux as computed per element and the flux projected onto the
   // smooth_flux_fes space.
   //double total_error = 0.0;
   errors_norm.SetSize(xfes->GetNE());
   for (int i = 0; i < xfes->GetNE(); i++)
   {
      errors_norm(i) = ComputeElementLpDistance(norm_p, i, smooth_flux, flux);
   }

   //double glob_error;
   //MPI_Allreduce(&total_error, &glob_error, 1, MPI_DOUBLE, MPI_SUM,
   //              xfes->GetComm());

   //return pow(glob_error, 1.0/norm_p);
   return fail;
}

