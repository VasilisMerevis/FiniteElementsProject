using System;

namespace FiniteElementsProject
{
    public class Element1D : IElements1D
    {
        protected double E, A;
        protected double[] nodesX, nodesY;
        protected double[,] localStiffnessMatrix, massMatrix;
        protected double[,] lambdaMatrix;
        public double[,] globalStiffnessMatrix { get; set; }
        protected double[] node1GlobalDisplacementVector, node2GlobalDisplacementVector;
        protected double[] node1XYInitial, node2XYInitial;
        protected double[] node1XYCurrent, node2XYCurrent;
        protected double cosInitial, sinInitial, lengthInitial;
        protected double cosCurrent, sinCurrent, lengthCurrent;
        public double[] internalGlobalForcesVector;

        public Element1D (double E, double A, double[] nodesX, double[] nodesY)
        {
            this.E = E;
            this.A = A;
            this.nodesX = nodesX;
            this.nodesY = nodesY;
            node1XYInitial = new[] { nodesX[0], nodesY[0] };
            node2XYInitial = new[] { nodesX[1], nodesY[1] };
            node1XYCurrent = new double[2];
            node2XYCurrent = new double[2];
        }

        #region Calculate_Geometric_Data
        protected void CalculateCurrentNodalCoordinates()
        {
            node1XYCurrent = new[]
            {
                node1XYInitial[0] + node1GlobalDisplacementVector[0],
                node1XYInitial[1] + node1GlobalDisplacementVector[1]
            };

            node2XYCurrent = new[]
            {
                node2XYInitial[0] + node2GlobalDisplacementVector[0],
                node2XYInitial[1] + node2GlobalDisplacementVector[1]
            };
        }

        protected double CalculateElementLength(double[] node1XY, double[] node2XY)
        {
            double length = Math.Sqrt(Math.Pow(node2XY[0] - node1XY[0], 2) + Math.Pow(node2XY[1] - node1XY[1], 2));
            return length;
        }

        protected double CalculateElementCosinus(double[] node1XY, double[] node2XY, double length)
        {
            double cosinus = (node2XY[0] - node1XY[0]) / length;
            return cosinus;
        }

        protected double CalculateElementSinus(double[] node1XY, double[] node2XY, double length)
        {
            double sinus = (node2XY[1] - node1XY[1]) / length;
            return sinus;
        }
        #endregion

        public virtual double[,] CreateLocalStiffnessMatrix()
        {
            return null;
        }

        public virtual double[,] CreateMassMatrix()
        {
            return null;
        }

        public virtual double[,] CreateLambdaMatrix()
        {
            return null;
        }

        public virtual double[,] CreateGlobalStiffnessMatrix()
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

