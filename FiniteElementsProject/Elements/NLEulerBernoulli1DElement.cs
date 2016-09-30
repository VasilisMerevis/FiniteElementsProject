using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class NLEulerBernoulli1DElement : Element1D
    {
        //double[] localDisplacementVector;
        //public double[] internalLocalForcesVector;
        private double I;
        //private double[,] Dmatrix;
        //private double[,] Bmatrix;
        private double betaAngleInitial;
        private double betaAngleCurrent;

        public NLEulerBernoulli1DElement(double E, double A, double I, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            this.I = I;
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
            this.internalGlobalForcesVector = new double[6];
        }

        private double CalculateElementBetaAngle(double[] node1XY, double[] node2XY)
        {
            double betaAngle = Math.Atan2( node2XY[1] - node1XY[1] , node2XY[0] - node1XY[0]);
            return betaAngle;
        }

        //private double[] CalculateLocalDisplacementVector()
        //{
        //    localDisplacementVector = new[]
        //    {
        //        lengthCurrent - lengthInitial,
        //        node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial,
        //        node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial
        //    };
        //    return localDisplacementVector;
        //}
      
        //private double[,] CreateDMatrix()
        //{
        //    Dmatrix = new[,]
        //    {
        //        { E * A / lengthInitial, 0, 0 },
        //        { 0 , 4 * E * I / lengthInitial, 2 * E * I / lengthInitial },
        //        { 0 , 2 * E * I / lengthInitial, 4 * E * I / lengthInitial }
        //    };

        //    return Dmatrix;
        //}

        //private double[,] CreateBMatrix()
        //{
        //    Bmatrix = new[,]
        //    {
        //        { -cosCurrent, -sinCurrent, 0, cosCurrent, sinCurrent, 0 },
        //        { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 1, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 0 },
        //        { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 0, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 1 }
        //    };

        //    return Bmatrix;
        //}

        //private double[] CalculateInternalLocalForcesVector()
        //{
        //    internalLocalForcesVector = VectorOperations.MatrixVectorProduct(Dmatrix, localDisplacementVector);
        //    return internalLocalForcesVector;
        //}

        private double[] CalculateInternalGlobalForcesVector()
        {
            //internalGlobalForcesVector = VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Bmatrix), internalLocalForcesVector);
            double N = A * E * (lengthCurrent - lengthInitial) / lengthInitial;
            double M1 = 2 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 4 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double M2 = 4 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 2 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double[] intforceV = new double[6];
            intforceV[0] = -(M2 * sinCurrent / lengthCurrent) - (M1 * sinCurrent / lengthCurrent) - N * cosCurrent;
            intforceV[1] = (M2 * cosCurrent / lengthCurrent) + (M1 * cosCurrent / lengthCurrent) - N * sinCurrent;
            intforceV[2] = M1;
            intforceV[3] = (M2 * sinCurrent / lengthCurrent) + (M1 * sinCurrent / lengthCurrent) + N * cosCurrent;
            intforceV[4] = -(M2 * cosCurrent / lengthCurrent) - (M1 * cosCurrent / lengthCurrent) + N * sinCurrent;
            intforceV[5] = M2;
            //internalGlobalForcesVector = intforceV;
            return intforceV; 
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double N = A * E * (lengthCurrent - lengthInitial) / lengthInitial;
            double M1 = 2 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 4 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double M2 = 4 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 2 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double sinb = sinCurrent;
            double cosb = cosCurrent;
            double cosb2 = cosCurrent * cosCurrent;
            double sinb2 = sinCurrent * sinCurrent;
            double cosbsinb = cosCurrent * sinCurrent;
            double Lc = lengthCurrent;
            double L0 = lengthInitial;
            double Lc2 = lengthCurrent * lengthCurrent;

            double[,] localStiffnessMatrix = new double[6, 6];
            localStiffnessMatrix[0, 0] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            localStiffnessMatrix[0, 1] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            localStiffnessMatrix[0, 2] = -6 * E * I * sinb / (L0 * Lc);
            localStiffnessMatrix[0, 3] = -N * sinb2 / Lc - 12 * E * I * sinb2 / (L0 * Lc2) + 2 * (M2 + M1) * cosbsinb / Lc2 - A * E * cosb2 / L0;
            localStiffnessMatrix[0, 4] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            localStiffnessMatrix[0, 5] = -6 * E * I * sinb / (L0 * Lc);

            localStiffnessMatrix[1, 0] = localStiffnessMatrix[0, 1];
            localStiffnessMatrix[1, 1] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            localStiffnessMatrix[1, 2] = 6 * E * I * cosb / (L0 * Lc);
            localStiffnessMatrix[1, 3] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            localStiffnessMatrix[1, 4] = -A * E * sinb2 / L0 - 2 * (M2 + M1) * cosbsinb / Lc2 - N * cosb2 / Lc - 12 * E * I * cosb2 / (L0 * Lc2);
            localStiffnessMatrix[1, 5] = 6 * E * I * cosb / (L0 * Lc);

            localStiffnessMatrix[2, 0] = localStiffnessMatrix[0, 2];
            localStiffnessMatrix[2, 1] = localStiffnessMatrix[1, 2];
            localStiffnessMatrix[2, 2] = 4 * E * I / L0;
            localStiffnessMatrix[2, 3] = 6 * E * I * sinb / (L0 * Lc);
            localStiffnessMatrix[2, 4] = -6 * E * I * cosb / (L0 * Lc);
            localStiffnessMatrix[2, 5] = 2 * E * I / L0;

            localStiffnessMatrix[3, 0] = localStiffnessMatrix[0, 3];
            localStiffnessMatrix[3, 1] = localStiffnessMatrix[1, 3];
            localStiffnessMatrix[3, 2] = localStiffnessMatrix[2, 3];
            localStiffnessMatrix[3, 3] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            localStiffnessMatrix[3, 4] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            localStiffnessMatrix[3, 5] = 6 * E * I * sinb / (L0 * Lc);

            localStiffnessMatrix[4, 0] = localStiffnessMatrix[0, 4];
            localStiffnessMatrix[4, 1] = localStiffnessMatrix[1, 4];
            localStiffnessMatrix[4, 2] = localStiffnessMatrix[2, 4];
            localStiffnessMatrix[4, 3] = localStiffnessMatrix[3, 4];
            localStiffnessMatrix[4, 4] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            localStiffnessMatrix[4, 5] = -6 * E * I * cosb / (L0 * Lc);

            localStiffnessMatrix[5, 0] = localStiffnessMatrix[0, 5];
            localStiffnessMatrix[5, 1] = localStiffnessMatrix[1, 5];
            localStiffnessMatrix[5, 2] = localStiffnessMatrix[2, 5];
            localStiffnessMatrix[5, 3] = localStiffnessMatrix[3, 5];
            localStiffnessMatrix[5, 4] = localStiffnessMatrix[4, 5];
            localStiffnessMatrix[5, 5] = 4 * E * I / L0;

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

            //Bmatrix = CreateBMatrix();
            //localDisplacementVector = CalculateLocalDisplacementVector();
            //Dmatrix = CreateDMatrix();
            //internalLocalForcesVector = CalculateInternalLocalForcesVector();
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
            //Bmatrix = CreateBMatrix();
            //localDisplacementVector = CalculateLocalDisplacementVector();
            //Dmatrix = CreateDMatrix();
            //internalLocalForcesVector = CalculateInternalLocalForcesVector();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            localStiffnessMatrix = CreateLocalStiffnessMatrix();
            globalStiffnessMatrix = localStiffnessMatrix;
        }
    }
}
