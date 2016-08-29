using System;

namespace FiniteElementsProject
{
	public class Element1D : IElements1D
	{
		protected double E, A;
		protected double[] nodesX, nodesY;
		protected double[,] localStiffnessMatrix;
		protected double[,] lambdaMatrix;
		public double[,] globalStiffnessMatrix { get; set; }
        protected double[] node1GlobalDisplacementVector, node2GlobalDisplacementVector;


        public Element1D (double E, double A, double[] nodesX, double[] nodesY)
		{
			this.E = E;
			this.A = A;
			this.nodesX = nodesX;
			this.nodesY = nodesY;
		}

		public virtual double[,] CreateLocalStiffnessMatrix()
		{
			return null;
		}

		public virtual double[,] CreateLambdaMatrix()
		{
			return null;
		}

		public double[,] CreateGlobalStiffnessMatrix()
		{
			double[,] lambdaTransposeMatrix = MatrixOperations.Transpose(lambdaMatrix);
			double[,] localStiffByLambda = MatrixOperations.MatrixProduct(localStiffnessMatrix, lambdaMatrix);
			globalStiffnessMatrix = MatrixOperations.MatrixProduct(lambdaTransposeMatrix, localStiffByLambda);

			return globalStiffnessMatrix;
		}

        public virtual void CalculateInitialValues()
        {
            CreateLocalStiffnessMatrix();
            CreateLambdaMatrix();
            CreateGlobalStiffnessMatrix();
        }

        public virtual void CalculateCurrentValues()
        {

        }

        public void SetGlobalDisplacementVector(double[] node1GlobalDisplacementVector, double[] node2GlobalDisplacementVector)
        {
            this.node1GlobalDisplacementVector = node1GlobalDisplacementVector;
            this.node2GlobalDisplacementVector = node2GlobalDisplacementVector;
        }
    }
}

