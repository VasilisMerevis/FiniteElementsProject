using System;

namespace FiniteElementsProject
{
	public class IterativeSolver : LinearSolver
	{
        private int maxIterations;
        private double tolerance = 1e-12;

        public IterativeSolver (double[,] stiffnessMatrix, double[] forceVector, int maxIterations)
		{
            this.stiffnessMatrix = stiffnessMatrix;
            this.forceVector = forceVector;
            solutionVector = new double[forceVector.Length];
			this.maxIterations = maxIterations;
		}

        override public void SolveWithMethod(string method)
        {
            switch (method)
            {
                case "PCG":
                    PCG();
                    break;
                case "Jacobi":
                    Jacobi();
                    break;
            }
        }

        private double [] PCG()
		{
			double [,] preconditioner = new double[stiffnessMatrix.GetLength(0), stiffnessMatrix.GetLength(1)];
			for (int i = 0; i < preconditioner.GetLength(0); i++) 
			{
				preconditioner [i, i] = 1 / stiffnessMatrix [i, i];
			}
			double[] residual = VectorOperations.VectorVectorSubtraction(
				forceVector,
				VectorOperations.MatrixVectorProduct (stiffnessMatrix, solutionVector)
				);
			double[] preconVector = VectorOperations.MatrixVectorProduct (preconditioner, residual);
			for (int iter = 0; iter < maxIterations; iter++) 
			{
				double[] u = VectorOperations.MatrixVectorProduct (stiffnessMatrix, preconVector);
				double residualDotOld = VectorOperations.VectorDotProduct (residual, preconVector);
				double alpha = residualDotOld / VectorOperations.VectorDotProduct (preconVector, u);
				solutionVector = VectorOperations.VectorVectorAddition(solutionVector,VectorOperations.VectorScalarProduct(preconVector, alpha));
				residual = VectorOperations.VectorVectorSubtraction (residual, VectorOperations.VectorScalarProduct (u, alpha));
				if (VectorOperations.VectorDotProduct (residual, residual) < tolerance) 
				{
					break;
				}
				double residualDotNew = VectorOperations.VectorDotProduct (residual, VectorOperations.MatrixVectorProduct(preconditioner,residual));
				double beta = residualDotNew / residualDotOld;
				preconVector = VectorOperations.VectorVectorAddition (
					VectorOperations.MatrixVectorProduct (preconditioner, residual),
					VectorOperations.VectorScalarProduct (preconVector, beta)
					);
			}
			return solutionVector;
		}

        private double[] Jacobi()
        {
            double[] oldSolution = new double[forceVector.Length];
            double[,] Diag = MatrixOperations.PutZerosInDiag(stiffnessMatrix);
            MatrixOperations.InvertDiagMatrix(Diag);
            double[,] zeroDiagStiff = MatrixOperations.GetDiagMatrix(stiffnessMatrix);
            for (int i = 0; i < maxIterations; i++)
            {
                oldSolution = solutionVector;
                solutionVector = VectorOperations.MatrixVectorProduct(
                    Diag,
                    VectorOperations.VectorVectorSubtraction(
                        forceVector,
                        VectorOperations.MatrixVectorProduct(zeroDiagStiff, solutionVector)
                        ));
                if (VectorOperations.VectorDotProduct(solutionVector, oldSolution) < tolerance)
                {
                    break;
                }
            }
            return solutionVector;
        }


    }
}

