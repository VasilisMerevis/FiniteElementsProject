using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class NLEulerBernoulli1DElement : Element1D
    {
        //int node1ID, node2ID;
        double[] node1XYInitial, node2XYInitial;
        double[] node1XYCurrent, node2XYCurrent;
        //private double[] node1GlobalDisplacementVector, node2GlobalDisplacementVector;
        double[] localDisplacementVector;
        public double[] internalLocalForcesVector, internalGlobalForcesVector;
        private double I;
        private double[,] Dmatrix;
        private double[,] Bmatrix;
        private double cosInitial, sinInitial, lengthInitial, betaAngleInitial;
        private double cosCurrent, sinCurrent, lengthCurrent, betaAngleCurrent;

        public NLEulerBernoulli1DElement(double E, double A, double I, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            this.I = I;
            //this.node1ID = localNode1;
            //this.node2ID = localNode2;
            node1XYInitial = new[] { nodesX[0], nodesY[0] };
            node2XYInitial = new[] { nodesX[1], nodesY[1] };
            node1XYCurrent = new double[2];
            node2XYCurrent = new double[2];
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
            this.internalGlobalForcesVector = new double[6];
        }

        //public void SetGlobalDisplacementVector(double[] node1GlobalDisplacementVector, double[] node2GlobalDisplacementVector)
        //{
        //    this.node1GlobalDisplacementVector = node1GlobalDisplacementVector;
        //    this.node2GlobalDisplacementVector = node2GlobalDisplacementVector;
        //}

        private void CalculateCurrentNodalCoordinates()
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

        private double CalculateElementLength(double[] node1XY, double[] node2XY)
        {
            double length = Math.Sqrt(Math.Pow(node2XY[0] - node1XY[0], 2) + Math.Pow(node2XY[1] - node1XY[1], 2));
            return length;
        }

        private double CalculateElementCosinus(double[] node1XY, double[] node2XY, double length)
        {
            double cosinus = (node2XY[0] - node1XY[0]) / length;
            return cosinus;
        }

        private double CalculateElementSinus(double[] node1XY, double[] node2XY, double length)
        {
            double sinus = (node2XY[1] - node1XY[1]) / length;
            return sinus;
        }

        private double CalculateElementBetaAngle(double[] node1XY, double[] node2XY)
        {
            double betaAngle = Math.Atan2( node2XY[1] - node1XY[1] , node2XY[0] - node1XY[0]);
            return betaAngle;
        }

        private double[] CalculateLocalDisplacementVector()
        {
            localDisplacementVector = new[]
            {
                lengthCurrent - lengthInitial,
                node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial,
                node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial
            };
            return localDisplacementVector;
        }
      

        private double[,] CreateDMatrix()
        {
            Dmatrix = new[,]
            {
                { E * A / lengthInitial, 0, 0 },
                { 0 , 4 * E * I / lengthInitial, 2 * E * I / lengthInitial },
                { 0 , 2 * E * I / lengthInitial, 4 * E * I / lengthInitial }
            };

            return Dmatrix;
        }

        private double[,] CreateBMatrix()
        {
            Bmatrix = new[,]
            {
                { -cosCurrent, -sinCurrent, 0, cosCurrent, sinCurrent, 0 },
                { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 1, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 0 },
                { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 0, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 1 }
            };

            return Bmatrix;
        }

        private double[] CalculateInternalLocalForcesVector()
        {
            internalLocalForcesVector = VectorOperations.MatrixVectorProduct(Dmatrix, localDisplacementVector);
            return internalLocalForcesVector;
        }

        private double[] CalculateInternalGlobalForcesVector()
        {
            internalGlobalForcesVector = VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Bmatrix), internalLocalForcesVector);
            return internalGlobalForcesVector; 
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double[,] transposeBmatrix = MatrixOperations.Transpose(Bmatrix);
            double[] zVector = new[] {sinCurrent, -cosCurrent, 0, -sinCurrent, cosCurrent, 0};
            double[] vVector = new[] {-cosCurrent, -sinCurrent, 0, cosCurrent, sinCurrent, 0};
            double[,] zzMatrix = VectorOperations.VectorVectorTensorProduct(zVector, zVector);
            double[,] vzMatrix = VectorOperations.VectorVectorTensorProduct(vVector, zVector);
            double[,] zvMatrix = VectorOperations.VectorVectorTensorProduct(zVector, vVector);
            double[,] vzPluszv = MatrixOperations.MatrixAddition(vzMatrix, zvMatrix);
            double[,] firstMember = MatrixOperations.ScalarMatrixProduct(internalLocalForcesVector[0] / lengthCurrent, zzMatrix);
            double[,] secondMember = MatrixOperations.ScalarMatrixProduct((internalLocalForcesVector[1] + internalLocalForcesVector[2]) / Math.Pow(lengthCurrent, 2), vzPluszv);
            double[,] thirdMember = MatrixOperations.MatrixProduct(transposeBmatrix, MatrixOperations.MatrixProduct(Dmatrix, Bmatrix));
            double[,] localStiffnessMatrix = MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(firstMember, secondMember), thirdMember);

            return localStiffnessMatrix;
        }

        public override double[,] CreateLambdaMatrix()
        {
            double length = Math.Sqrt(Math.Pow((nodesY[1] - nodesY[0]), 2) + Math.Pow((nodesX[1] - nodesX[0]), 2));
            double cosine = (nodesX[1] - nodesX[0]) / length;
            double sine = (nodesY[1] - nodesY[0]) / length;
            lambdaMatrix = new[,]
            { { cosine, sine, 0, 0, 0, 0 },
            { -sine, cosine, 0, 0, 0, 0 },
            { 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, cosine, sine, 0 },
            { 0, 0, 0, -sine, cosine, 0 },
            { 0, 0, 0, 0, 0, 1 }};

            return lambdaMatrix;
        }

        public override void CalculateInitialValues()
        {
            lengthInitial = CalculateElementLength(node1XYInitial, node2XYInitial);
            sinInitial = CalculateElementSinus(node1XYInitial, node2XYInitial, lengthInitial);
            cosInitial = CalculateElementCosinus(node1XYInitial, node2XYInitial, lengthInitial);
            betaAngleInitial = CalculateElementBetaAngle(node1XYInitial, node2XYInitial);

            lengthCurrent = lengthInitial;
            sinCurrent = sinInitial;
            cosCurrent = cosInitial;
            betaAngleCurrent = betaAngleInitial;

            Bmatrix = CreateBMatrix();
            localDisplacementVector = CalculateLocalDisplacementVector();
            Dmatrix = CreateDMatrix();
            internalLocalForcesVector = CalculateInternalLocalForcesVector();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            localStiffnessMatrix = CreateLocalStiffnessMatrix();
            globalStiffnessMatrix = localStiffnessMatrix;

        }

        public override void CalculateCurrentValues()
        {
            CalculateCurrentNodalCoordinates();
            lengthCurrent = CalculateElementLength(node1XYCurrent, node2XYCurrent);
            sinCurrent = CalculateElementSinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
            cosCurrent = CalculateElementCosinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
            betaAngleCurrent = CalculateElementBetaAngle(node1XYCurrent, node2XYCurrent);
            Bmatrix = CreateBMatrix();
            localDisplacementVector = CalculateLocalDisplacementVector();
            Dmatrix = CreateDMatrix();
            internalLocalForcesVector = CalculateInternalLocalForcesVector();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            localStiffnessMatrix = CreateLocalStiffnessMatrix();
            globalStiffnessMatrix = localStiffnessMatrix;
        }
    }
}
