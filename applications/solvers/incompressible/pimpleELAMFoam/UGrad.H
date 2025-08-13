//----------------------------------------------------------------------------------------------------------------------------------
// Defining the velocity gradient tensor/mean velocity gradient tensor as a new variable to access below


const volVectorField& USource = mesh.foundObject<volVectorField>("UMean") ? mesh.lookupObject<volVectorField>("UMean") : mesh.lookupObject<volVectorField>("U");


//volTensorField gradU = fvc::grad(UMean);
volTensorField gradU = fvc::grad(USource);

// Creating variables for  each component of the velocity gradient tensor
volScalarField dUdx=gradU.component(0);
volScalarField dVdx=gradU.component(1);
volScalarField dWdx=gradU.component(2);
volScalarField dUdy=gradU.component(3);
volScalarField dVdy=gradU.component(4);
volScalarField dWdy=gradU.component(5);
volScalarField dUdz=gradU.component(6);
volScalarField dVdz=gradU.component(7);
volScalarField dWdz=gradU.component(8);

// Creating variables for each component of the velocity vector 
volScalarField Ux = USource.component(0);
volScalarField Uy = USource.component(1);
volScalarField Uz = USource.component(2);

// Looping through each cell volume
forAll (U, celli)
{

	NodeNo[celli] = celli;

	Ax[celli] = ( dUdx[celli] * Ux[celli] ) + ( dUdy[celli] * Uy[celli] ) + ( dUdz[celli] * Uz[celli] );
	Ay[celli] = ( dVdx[celli] * Ux[celli] ) + ( dVdy[celli] * Uy[celli] ) + ( dVdz[celli] * Uz[celli] );
	Az[celli] = ( dWdx[celli] * Ux[celli] ) + ( dWdy[celli] * Uy[celli] ) + ( dWdz[celli] * Uz[celli] );
	AcclM[celli] = Foam::pow(((Ax[celli] * Ax[celli]) + (Ay[celli] * Ay[celli]) + (Az[celli] * Az[celli])),0.5);
	STRXYZUVW[celli] = abs(dUdx[celli]) + abs(dUdy[celli]) + abs(dUdz[celli]) + abs(dVdx[celli]) + abs(dVdy[celli]) \
			   	+ abs(dVdz[celli]) + abs(dWdx[celli]) + abs(dWdy[celli]) + abs(dWdz[celli]);
}

// Updating boundary conditions
NodeBC.correctBoundaryConditions();
//WSC.correctBoundaryConditions();
Ax.correctBoundaryConditions();
Ay.correctBoundaryConditions();
Az.correctBoundaryConditions();
AcclM.correctBoundaryConditions();
STRXYZUVW.correctBoundaryConditions();
NodeNo.correctBoundaryConditions();
//-------------------------------------------------------------------------------------------------------------------------------------


/*
// Defining the velocity gradient tensor/mean velocity gradient tensor as a new variable to access below
volTensorField gradU = fvc::grad(U);
volTensorField gradUMean = fvc::grad(UMean);

// Creating variables for  each component of the velocity gradient tensor
volScalarField dUdx=gradU.component(0);
volScalarField dVdx=gradU.component(1);
volScalarField dWdx=gradU.component(2);
volScalarField dUdy=gradU.component(3);
volScalarField dVdy=gradU.component(4);
volScalarField dWdy=gradU.component(5);
volScalarField dUdz=gradU.component(6);
volScalarField dVdz=gradU.component(7);
volScalarField dWdz=gradU.component(8);

// Creating variables for  each component of the mean velocity gradient tensor
volScalarField dUMdx=gradUMean.component(0);
volScalarField dVMdx=gradUMean.component(1);
volScalarField dWMdx=gradUMean.component(2);
volScalarField dUMdy=gradUMean.component(3);
volScalarField dVMdy=gradUMean.component(4);
volScalarField dWMdy=gradUMean.component(5);
volScalarField dUMdz=gradUMean.component(6);
volScalarField dVMdz=gradUMean.component(7);
volScalarField dWMdz=gradUMean.component(8);


// Creating variables for each component of the velocity vector 
volScalarField Ux = U.component(0);
volScalarField Uy = U.component(1);
volScalarField Uz = U.component(2);

// Creating variables for each component of the mean velocity vector 
volScalarField UMeanx = UMean.component(0);
volScalarField UMeany = UMean.component(1);
volScalarField UMeanz = UMean.component(2);


// Looping through each cell volume
forAll (alpha1,celli)
{

	if (mesh.foundObject<volVectorField>("UMean")) // check to see if UMean exists in the mesh & calculates accelerations if so
	{
	Ax[celli] = ( dUMdx[celli] * UMeanx[celli] ) + ( dUMdy[celli] * UMeany[celli] ) + ( dUMdz[celli] * UMeanz[celli] );
	Ay[celli] = ( dVMdx[celli] * UMeanx[celli] ) + ( dVMdy[celli] * UMeany[celli] ) + ( dVMdz[celli] * UMeanz[celli] );
	Az[celli] = ( dWMdx[celli] * UMeanx[celli] ) + ( dWMdy[celli] * UMeany[celli] ) + ( dWMdz[celli] * UMeanz[celli] );
	
	AcclM[celli] = Foam::pow(((Ax[celli] * Ax[celli]) + (Ay[celli] * Ay[celli]) + (Az[celli] * Az[celli])),0.5);
	STRXYZUVW[celli] = abs(dUdx[celli]) + abs(dUdy[celli]) + abs(dUdz[celli]) + abs(dVdx[celli]) + abs(dVdy[celli]) \
			   + abs(dVdz[celli]) + abs(dWdx[celli]) + abs(dWdy[celli]) + abs(dWdz[celli]);
	}
	else // caclulates accelerations using standard velocity vectors if UMean does not exist
	{	
	Ax[celli] = ( dUdx[celli] * Ux[celli] ) + ( dUdy[celli] * Uy[celli] ) + ( dUdz[celli] * Uz[celli] );
	Ay[celli] = ( dVdx[celli] * Ux[celli] ) + ( dVdy[celli] * Uy[celli] ) + ( dVdz[celli] * Uz[celli] );
	Az[celli] = ( dWdx[celli] * Ux[celli] ) + ( dWdy[celli] * Uy[celli] ) + ( dWdz[celli] * Uz[celli] );
	
	AcclM[celli] = Foam::pow(((Ax[celli] * Ax[celli]) + (Ay[celli] * Ay[celli]) + (Az[celli] * Az[celli])),0.5);
	STRXYZUVW[celli] = abs(dUdx[celli]) + abs(dUdy[celli]) + abs(dUdz[celli]) + abs(dVdx[celli]) + abs(dVdy[celli]) \
			   + abs(dVdz[celli]) + abs(dWdx[celli]) + abs(dWdy[celli]) + abs(dWdz[celli]);
	}
	

	if (alpha1[celli]>0.5)
	{

		NodeBC[celli]=0;
		WSC[celli]=-1;

	}	
}

// Updating boundary conditions
NodeBC.correctBoundaryConditions();
WSC.correctBoundaryConditions();
Ax.correctBoundaryConditions();
Ay.correctBoundaryConditions();
Az.correctBoundaryConditions();
AcclM.correctBoundaryConditions();
STRXYZUVW.correctBoundaryConditions();
*/
