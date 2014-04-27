/* Bevan, Josh UMass Lowell 2013-2014
 *Developed with the help of Prof. Trelles, UMass Lowell
 *In partial satisfaction of directed study 22.602 Fall 2013
 *Solves the 2D Laplacian on an unstructured grid with KSP.*/
static char help[] = "Solves the 2D Laplacian on an unstructured grid with KSP. Bevan 2014\n\n";

#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x, b;             /* approx solution, RHS*/
  Mat            A;                /* linear system matrix */
  KSP            ksp;              /* linear solver context */
  PC             pc;               /* preconditioner context */
  PetscReal      norm; 			   /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       n = 81;				/*Number of DOFs*/
  /*For simplicity of use the mesh input file lengths need to be manually input here JJBHot[],JJBCold[],JJBElems[], and JJBNodes[]*/
  PetscInt		   row[3],col[3],its,elnd[3],JJBElems[384],iter=0,NbElems,NbVertices,numint,xx=0,yy=1,ie,rstart,rend,nlocal,NbHot,NbCold;
  PetscInt 	   JJBHot[9],JJBCold[9],start,end;
  PetscScalar    value[9],n1[2],n2[2],n3[2],gn1[2],gn2[2],gn3[2],A11,A22,A33,A12,A13,A23,JJBNodes[162],numscal, TriArea;
  PetscScalar	   HotBC[9], ColdBC[9];	/*Dirichlet BC values, top and bottom domains implicitly Neumann=0*/
  PetscViewer	   viewer;

  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = PetscOptionsGetInt(NULL,"-sizer",&n,NULL);CHKERRQ(ierr);
/*---------------------------------------------------------------------------------------*/
  /* Create vectors. Create one vector and duplicate as needed.
   * The second argument to VecSetSizes() below causes PETSc to decide
   * how many elements per processor are assigned */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  /* Find start and end mesh points on each processor for the 
   * interior of the mesh. This partitioning is based upon
   * the choices made for the vector in VecSetSizes()*/
  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);
/*---------------------------------------------------------------------------------------*/
  /* Create matrix A with MatCreate(); matrix format can be 
   * specified at runtime
   * Use nlocal as the local size of the matrix, this ensures it will
   * align with the vectors from above */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
/*---------------------------------------------------------------------------------------*/
/*Read in mesh */
  
	/*Read in elem/node associations*/
	FILE *file = fopen("JJBElems", "r");
	while(fscanf(file, "%d", &numint) == 1) {
		JJBElems[iter] = numint;
		iter++;
	}
	NbElems = iter/3;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: Read in %D elements.\n",NbElems);CHKERRQ(ierr);
	fclose(file);

	/*Read in vertex coordinates*/
	FILE *file2 = fopen("JJBNodes", "r");
	iter=0;
	while(fscanf(file2, "%lf", &numscal) == 1) { /*Note %lf for doubles to match PetscScalar */
		JJBNodes[iter] = numscal;
		iter++;
	}
	NbVertices = iter/2;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: Read in %D nodes.\n",NbVertices);CHKERRQ(ierr);
	fclose(file2);
	
	/*Read in Hot BC*/
	FILE *file3 = fopen("JJBHot", "r");
	iter=0;
	while(fscanf(file3, "%d", &numint) == 1) {
		JJBHot[iter] = numint-1;
		iter++;
	}
	NbHot = iter;
	fclose(file3);
	
	/*Read in Cold BC*/
	FILE *file4 = fopen("JJBCold", "r");
	iter=0;
	while(fscanf(file4, "%d", &numint) == 1) {
		JJBCold[iter] = numint-1;
		iter++;
	}
	NbCold = iter;
	fclose(file4);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: Read in %D explicit boundary nodes.\n",NbHot+NbCold);CHKERRQ(ierr);

/*---------------------------------------------------------------------------------------*/
  /* Assemble matrix.

     The linear system is distributed across the processors by
     chunks of contiguous rows, which correspond to contiguous
     sections of the mesh on which the problem is discretized.
     For matrix assembly, each processor contributes entries for
     the part that it owns locally.  */
	 
  /* From manual: The routine MatSetValuesBlocked() may offer much
   * better efficiency for users of block sparse formats
   * (MATSEQBAIJ and MATMPIBAIJ) */
	ierr = MatGetOwnershipRange(A, &start, &end);CHKERRQ(ierr);
	if (end!=n){end -= 1;}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: ____Processor Partitioning____\n");CHKERRQ(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"JJB: (Nodes) Start: %D, End %D\n",start,end);CHKERRQ(ierr);
	if (start!=0){start = (floor((((1.0*start)-1)/n)*NbElems)+1)*3;}
	end = floor(((1.0*end)/n)*NbElems)*3;
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"JJB: (Elements) Start: %D, End %D\n",start,end);CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

for (ie=start;(ie<end+1)&&(ie<NbElems*3);ie+=3){
	/*For each element find the associated node numbers */
    elnd[0] = JJBElems[ie]-1;
    elnd[1] = JJBElems[ie+1]-1;
    elnd[2] = JJBElems[ie+2]-1;
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"JJB: ie %D\n",ie);CHKERRQ(ierr);
    
	/* Determine element node locations */
	n1[xx] = JJBNodes[2*elnd[0]];
	n1[yy] = JJBNodes[(2*elnd[0])+1];
	n2[xx] = JJBNodes[2*elnd[1]];
	n2[yy] = JJBNodes[(2*elnd[1])+1];
	n3[xx] = JJBNodes[2*elnd[2]];
	n3[yy] = JJBNodes[(2*elnd[2])+1];
	
    /*Calculate bounded area*/
    TriArea = fabs(( (n1[xx]*(n2[yy]-n3[yy])) + (n2[xx]*(n3[yy]-n1[yy])) + (n3[xx]*(n1[yy]-n2[yy])) )/2);

    /*Calculate gradients:
	 * without the benefit of a low overhead linear solver these
	 * were algebraicly solved. The equivalent Matlab code is:
	 Plane = [X0 Y0 1;...
             X1 Y1 1;...
             X2 Y2 1];
    Gradient(:,0) = Plane\[1;0;0];
    Gradient(:,1) = Plane\[0;1;0]; 
    Gradient(:,2) = Plane\[0;0;1];
    Gradient(3,:) = 0; %Remove unwanted c coeff */
	 
	if (n2[yy]!=n3[yy]){
		gn1[xx] = -1/( ( (n2[xx]-n3[xx]) * (n1[yy]-n3[yy]) / (n2[yy]-n3[yy]) ) - (n1[xx]-n3[xx]) );
	}else{
		gn1[xx] = 0;}
	if (n2[xx]!=n3[xx]){
		gn1[yy] = -1/( ( (n2[yy]-n3[yy]) * (n1[xx]-n3[xx]) / (n2[xx]-n3[xx]) ) - (n1[yy]-n3[yy]) );
	}else{
		gn1[yy] = 0;}

	if (n1[yy]!=n3[yy]){
		gn2[xx] = -1/( ( (n1[xx]-n3[xx]) * (n2[yy]-n3[yy]) / (n1[yy]-n3[yy]) ) - (n2[xx]-n3[xx]) );
	}else{
		gn2[xx] = 0;}
	if (n1[xx]!=n3[xx]){
		gn2[yy] = -1/( ( (n1[yy]-n3[yy]) * (n2[xx]-n3[xx]) / (n1[xx]-n3[xx]) ) - (n2[yy]-n3[yy]) );
	}else{
		gn2[yy] = 0;}

	if (n1[yy]!=n2[yy]){
		gn3[xx] = -1/( ( (n1[xx]-n2[xx]) * (n3[yy]-n2[yy]) / (n1[yy]-n2[yy]) ) - (n3[xx]-n2[xx]) );
	}else{
		gn3[xx] = 0;}
	if (n1[xx]!=n2[xx]){
		gn3[yy] = -1/( ( (n1[yy]-n2[yy]) * (n3[xx]-n2[xx]) / (n1[xx]-n2[xx]) ) - (n3[yy]-n2[yy]) );
	}else{
		gn3[yy] = 0;}
	
	/*Calculate local elemental matrix*/
	A11 = TriArea * (gn1[xx]*gn1[xx] + gn1[yy]*gn1[yy]);
	A22 = TriArea * (gn2[xx]*gn2[xx] + gn2[yy]*gn2[yy]);
	A33 = TriArea * (gn3[xx]*gn3[xx] + gn3[yy]*gn3[yy]);
	A12 = TriArea * (gn1[xx]*gn2[xx] + gn1[yy]*gn2[yy]);
	A13 = TriArea * (gn1[xx]*gn3[xx] + gn1[yy]*gn3[yy]);
	A23 = TriArea * (gn2[xx]*gn3[xx] + gn2[yy]*gn3[yy]);
	
	/*Create value array to be stamped*/
	value[0] = A11;
	value[1] = A12;
	value[2] = A13;
	
	value[3] = A12;
	value[4] = A22;
	value[5] = A23;
	
	value[6] = A13;
	value[7] = A23;
	value[8] = A33;
	
	/*ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: ie: %D\n",ie);CHKERRQ(ierr);
	for (iter=0;iter<9;iter+=1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: Val: %lf\n",value[iter]);CHKERRQ(ierr);
	}*/
	
	/*Create global node index arrays for stamping*/
	row[0] = elnd[0];
	row[1] = elnd[1];
	row[2] = elnd[2];

	/*-*/
	col[0] = elnd[0];
	col[1] = elnd[1];
	col[2] = elnd[2];

	/* Stamp in local elemental matrix into global matrix of form: */ 
	/*					   
	A(N1, N1) = A(N1, N1) + A_elemental(1,1);
    A(N2, N1) = A(N2, N1) + A_elemental(2,1);
    A(N3, N1) = A(N3, N1) + A_elemental(3,1);
 
    A(N1, N2) = A(N1, N2) + A_elemental(1,2);
    A(N2, N2) = A(N2, N2) + A_elemental(2,2);
    A(N3, N2) = A(N3, N2) + A_elemental(3,2);
    
    A(N1, N3) = A(N1, N3) + A_elemental(1,3);
    A(N2, N3) = A(N2, N3) + A_elemental(2,3);
    A(N3, N3) = A(N3, N3) + A_elemental(3,3);*/
    ierr = MatSetValues(A,3,row,3,col,value,ADD_VALUES);CHKERRQ(ierr);
}

  /* Assemble the matrix */
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	/* Apply BCs
	* MatZeroRows allows easy application of Dirichlet BCs by zeroing all entries
	* except the main diag (i.e. the self-reference needed to enforce BCs */
	ierr = VecSet(b,0);CHKERRQ(ierr);

	for(iter=0;iter<NbHot;iter++){
		HotBC[iter] = 100;
	}
	for(iter=0;iter<NbCold;iter++){
		ColdBC[iter] = 10;
	}
	ierr   = MatZeroRows(A,NbHot,JJBHot,1,0,0);CHKERRQ(ierr);
	ierr = VecSetValues(b,NbHot,JJBHot,HotBC,INSERT_VALUES);CHKERRQ(ierr);

	ierr   = MatZeroRows(A,NbCold,JJBCold,1,0,0);CHKERRQ(ierr);
	ierr = VecSetValues(b,NbCold,JJBCold,ColdBC,INSERT_VALUES);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
	
	ierr   = MatView(A,0);CHKERRQ(ierr);

/*---------------------------------------------------------------------------------------*/
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr); /* Create linear solver context */

  /* Set operators. A also serves as preconditioning matrix */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /* Set linear solver defaults for problem.
     Extract KSP and PC contexts from the KSP context,
	 to directly call any KSP and PC routines to set options */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  /* Set runtime options,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr); 							/* Solve linear system */
/*---------------------------------------------------------------------------------------*/
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);  /* View solver info */

  /* Check the error*/
  ierr = KSPGetResidualNorm(ksp,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"JJB: Norm of error %G, Iterations %D\n",norm,its);CHKERRQ(ierr);
  
  /* Output solution vector to file for external plotting*/
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"JJBx.output",&viewer);CHKERRQ(ierr);
  ierr = VecView(x,viewer);CHKERRQ(ierr);

  /* Free work space */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  /* PetscFinalize() before exiting program. Provides summary
   * and diagnostic information if certain runtime options 
   * are chosen (e.g., -log_summary) */
  ierr = PetscFinalize();
  return 0;
}