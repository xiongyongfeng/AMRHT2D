#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h" 
#include <vector>
//#include "petsc.h"
//#include <petscksp.h>
//#include "petscsys.h"
//#include "petscdef.h"
//#include <petscvec.h>
//#include <petscmat.h>
//#include <petscpc.h>
//#include <petscviewer.h>
//#include <petscis.h>


#define PI 3.1415926535897

using namespace std;

double f (double x, double y) { 
  return (cos(2.0 * PI * x) * sin(2.0 * PI * y));
}

double df (double x, double y) {
  return (-8.0 * PI * PI * f(x,y));
}

int main (int argc, char **args){
  
  cell * c;
  //int level = 0;
  int number_of_levels = 1;
  int nxb = 4;
  int nyb = 4;
  vector<int> * max_dimension_by_level;
  vector<double> *rhs, *u, * ue, *ff, *A;
  vector<int> *JA, *IA;

  /* PETSc 
  Vec x, b, v, erro; // approx solution, RHS, exact solution 
  Mat AP;       // linear system matrix 
  KSP ksp;     // linear solver context 
  PC pc;       // preconditioner context 
  PetscReal norminf, norm2; // norm of solution error 
  PetscInt ip, np, iter;
  PetscErrorCode ierr;
  PetscMPIInt size;  */
  
  //double delta_x, delta_y;
  dominio * D;
  int ncell = 0;
  int ncellv = 0;
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.;
  xend = yend = 1.0;
  D = new dominio (xbegin, ybegin, xend, yend);
  mesh * M;

  //delta_x = delta_y = 0.05;
  max_dimension_by_level = new vector<int>;
  rhs = new vector<double>;
  u = new vector<double>;
  ff = new vector<double>;
    
  //srand(time(NULL));
  srand(12345);
  
  /**********initiallize mesh*******************/
  for (int i = 0; i < number_of_levels; i++)
    max_dimension_by_level->push_back((nxb * pow(2,i)) * (nyb * pow(2,i)));
  M = new mesh(D, number_of_levels, nxb ,nyb, max_dimension_by_level);
  /*********************************************/
  /* Nivel Base - Nivel 0 */
  for (int y = 0; y < nxb; y++)
    for (int x = 0; x < nyb; x++){
      c = new cell(x, y, 0,-1);
      M->insert(c);
    }
    
  /**************Refinement*************************/
  list <cell *> * l = M->get_list_cell_by_level (0);
  //list <cell *>::iterator it = l->begin();
  //hash_table * Hnew = M->get_hash_table();

  /* Caso teste com refinamento estatico */
  /*while (it != l->end()) {
    if (((*it)->get_cell_x() >= 1 && (*it)->get_cell_x() < 2) && ((*it)->get_cell_y () >= 1 && (*it)->get_cell_y() <=2)) {
      c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
      assert (c != NULL);
      it = M->split_and_insert(c);
    }
    else
      it++;
      }*/

  /******************************************/
  /* Apos a criacao da malha gerar os coeficientes da matriz uma lista de pesos para cada celula da malha */

  M->get_hash_table()->print_information();
  
  M->create_unstructured_mesh(&f, &df);

  ncell = M->counting_mesh_cell();
  ncellv = M->neighbours_all_cell();
  cout << ncell << " " << ncellv << endl;
  vector <cell *> * elem;
  elem = new vector <cell *> (ncell);
  IA = new vector <int> (ncell+1,-1);
  JA = new vector <int> (ncellv,-1);
  A = new vector <double> (ncellv,0.0);
  u = new vector <double> (ncell,0.0);
  
  M->mesh_adress(ncell, elem);
  M->create_matrix_df (A, JA, IA, ncell, ncellv, elem);
  M->print_matrix(IA, JA, A, ncell, ncellv);
  
  //M->calculation_function (&df, ff);
  //M->calculation_function (&f, ue);
  M->rhs_dirichlet_boundary_conditions (&df, rhs);
  cout << "rhsmax = " << M->normmax(rhs, ncell) << endl;

  delete IA;
  delete A;
  delete JA;
  delete rhs;
  delete u;
  rhs = new vector <double> (5);
  u = new vector <double> (5, 0.0);
  ue = new vector <double> (5, 0.0);
  IA = new vector <int> (6,-1);
  JA = new vector <int> (12,-1);
  A = new vector <double> (12,0.0);
  for(int ik = 0; ik < 12; ik++)
    A->push_back(ik+1);
  JA->at(0) = 1; JA->at(1) = 4;
  JA->at(2) = 1; JA->at(3) = 2; JA->at(4) = 4;
  JA->at(5) = 1; JA->at(6) = 3; JA->at(7) = 4; JA->at(8) = 5;
  JA->at(9) = 3; JA->at(10) = 4;
  JA->at(11) = 5;

  IA->at(0) = 1; IA->at(1) = 3; IA->at(2) = 6; IA->at(3) = 10;
  IA->at(4) = 12; IA->at(5) = 13;
  ue->at(0) = 1; ue->at(1) = -1; ue->at(3) = -1;
  rhs->at(0) = -1; rhs->at(1) = -6; rhs->at(2) = -2; rhs->at(3) = -11;
  rhs->at(4) = 0;
  cout << "P" << endl;
  M->gs_csr(u, A, IA, JA, rhs, 5);
  cout << "|u - ue|_2 = " << M->erron2(u, ue, 5) << endl;
  
  /* PETSc Solver */
  //PetscInitialize(&argc, &args, (char *)0, PETSC_NULL);
  
  //ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  //if(size != 1)
  //  cout << "This is a uniprocessor example only!"<< endl;
  //ierr = PetscOptionsGetInt(NULL,NULL,"-n",&np,NULL);CHKERRQ(ierr);

  /* Matrix and right-hand-size: Ax = b 

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,ncell);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x, &b);CHKERRQ(ierr);
  ierr = VecDuplicate(x, &v);CHKERRQ(ierr);
  ierr = VecDuplicate(x, &erro);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&AP);CHKERRQ(ierr);
  ierr = MatSetSizes(AP,PETSC_DECIDE,PETSC_DECIDE,ncell,ncell);CHKERRQ(ierr);
  ierr = MatSetType(AP, MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(AP, 0, ncellv); CHKERRQ(ierr);
 
  ierr = MatSetFromOptions(AP);CHKERRQ(ierr);
  ierr = MatSetUp(AP);CHKERRQ(ierr);
 
  for(ip = 0; ip < ncell; ip++){
    ierr = VecSetValue(b, ip, rhs[ip], INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v, ip, ff[ip], INSERT_VALUES); CHKERRQ(ierr);
    for(jp = IA[ip]; jp < IA[ip+1]; jp++){
      ierr = MatSetValue(AP, ip, JA[jp], A[jp], INSERT_VALUES);CHKERRQ(ierr);
      cout << "i = " << ip << " j = " << JA[ip] << " Aij = " << A[jp] << endl;
    }
  }

  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(erro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(erro);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(AP, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(AP, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  /* Solver 
  
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, AP, AP);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPBCGS);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, 1.0e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /*Solution and data view 
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  KSPGetIterationNumber(ksp,iter);
 
  cout << "iteration = " << iter << endl;
  ierr = VecCopy(x, erro);CHKERRQ(ierr);
  ierr = VecAXPY(erro, -1.0, v);CHKERRQ(ierr);
  ierr = VecNorm(erro, NORM_INFINITY, norminf);CHKERRQ(ierr);
  ierr = VecNorm(erro, NORM_2, norm2);chkerrq(ierr);
  cout << "Norm_inf/Norm_2 = " << norminf << " " << norm2 << endl;
    
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&v); CHKERRQ(ierr);
  
  ierr = MatDestroy(&AP);CHKERRQ(ierr);*/

  //PetscFinalize();
  delete rhs;
  delete u;
  delete ff;
  delete elem;
  delete IA;
  delete A;
  delete JA;
	  
  return 0;
}
