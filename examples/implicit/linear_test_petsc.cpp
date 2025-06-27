#include <petscksp.h>

int main(int argc, char **args) {
    PetscInitialize(&argc, &args, nullptr, nullptr);

    Mat A;
    Vec b, x;
    KSP ksp;

    // Matrix size
    PetscInt n = 3;

    // Create matrix and vectors
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(A);
    MatSetUp(A);

    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, PETSC_DECIDE, n);
    VecSetFromOptions(b);

    VecDuplicate(b, &x); // x will hold the solution

    // Fill matrix A
    MatSetValue(A, 0, 0, 3.0, INSERT_VALUES);
    MatSetValue(A, 0, 1, 2.0, INSERT_VALUES);
    MatSetValue(A, 0, 2, -1.0, INSERT_VALUES);

    MatSetValue(A, 1, 0, 2.0, INSERT_VALUES);
    MatSetValue(A, 1, 1, -2.0, INSERT_VALUES);
    MatSetValue(A, 1, 2, 4.0, INSERT_VALUES);

    MatSetValue(A, 2, 0, -1.0, INSERT_VALUES);
    MatSetValue(A, 2, 1, 0.5, INSERT_VALUES);
    MatSetValue(A, 2, 2, -1.0, INSERT_VALUES);

    // Fill vector b
    VecSetValue(b, 0, 1.0, INSERT_VALUES);
    VecSetValue(b, 1, -2.0, INSERT_VALUES);
    VecSetValue(b, 2, 0.0, INSERT_VALUES);

    // Assemble matrix and vector
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // Create linear solver
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);

    // Solve Ax = b
    KSPSolve(ksp, b, x);

    // View the result
    PetscPrintf(PETSC_COMM_WORLD, "Solution x:\n");
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    // Clean up
    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

    PetscFinalize();
    return 0;
}
