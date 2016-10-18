using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class NLEulerBernoulli1DElement : Element1D
    {
        private readonly double I;
        private double betaAngleInitial;
        private double betaAngleCurrent;

        public NLEulerBernoulli1DElement(double E, double A, double I, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            this.I = I;
            node1GlobalDisplacementVector = new double[3];
            node2GlobalDisplacementVector = new double[3];
            globalStiffnessMatrix = new double[6, 6];
            internalGlobalForcesVector = new double[6];
        }

        private double CalculateElementBetaAngle(double[] node1XY, double[] node2XY)
        {
            double betaAngle = Math.Atan2( node2XY[1] - node1XY[1] , node2XY[0] - node1XY[0]);
            return betaAngle;
        }

        #region Local_Values_Calculations
        private double[] CalculateLocalDisplacementVector()
        {
            double[] localDisplacementVector = new[]
            {
                lengthCurrent - lengthInitial,
                node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial,
                node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial
            };
            return localDisplacementVector;
        }

        private double[] CalculateInternalLocalForcesVector()
        {
            double[,] Dmatrix = new[,]
            {
                { E * A / lengthInitial, 0, 0 },
                { 0 , 4 * E * I / lengthInitial, 2 * E * I / lengthInitial },
                { 0 , 2 * E * I / lengthInitial, 4 * E * I / lengthInitial }
            };

            double[] localDisplacementVector = CalculateLocalDisplacementVector();
            double[] internalLocalForcesVector = VectorOperations.MatrixVectorProduct(Dmatrix, localDisplacementVector);
            return internalLocalForcesVector;
        }
        #endregion

        #region Global_Values_Calculations
        private double[] CalculateInternalGlobalForcesVector()
        {
            
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
            
            return intforceV; 
        }

        public override double[,] CreateGlobalStiffnessMatrix()
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

            double[,] globalStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix[0, 0] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            globalStiffnessMatrix[0, 1] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            globalStiffnessMatrix[0, 2] = -6 * E * I * sinb / (L0 * Lc);
            globalStiffnessMatrix[0, 3] = -N * sinb2 / Lc - 12 * E * I * sinb2 / (L0 * Lc2) + 2 * (M2 + M1) * cosbsinb / Lc2 - A * E * cosb2 / L0;
            globalStiffnessMatrix[0, 4] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            globalStiffnessMatrix[0, 5] = -6 * E * I * sinb / (L0 * Lc);

            globalStiffnessMatrix[1, 0] = globalStiffnessMatrix[0, 1];
            globalStiffnessMatrix[1, 1] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[1, 2] = 6 * E * I * cosb / (L0 * Lc);
            globalStiffnessMatrix[1, 3] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            globalStiffnessMatrix[1, 4] = -A * E * sinb2 / L0 - 2 * (M2 + M1) * cosbsinb / Lc2 - N * cosb2 / Lc - 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[1, 5] = 6 * E * I * cosb / (L0 * Lc);

            globalStiffnessMatrix[2, 0] = globalStiffnessMatrix[0, 2];
            globalStiffnessMatrix[2, 1] = globalStiffnessMatrix[1, 2];
            globalStiffnessMatrix[2, 2] = 4 * E * I / L0;
            globalStiffnessMatrix[2, 3] = 6 * E * I * sinb / (L0 * Lc);
            globalStiffnessMatrix[2, 4] = -6 * E * I * cosb / (L0 * Lc);
            globalStiffnessMatrix[2, 5] = 2 * E * I / L0;

            globalStiffnessMatrix[3, 0] = globalStiffnessMatrix[0, 3];
            globalStiffnessMatrix[3, 1] = globalStiffnessMatrix[1, 3];
            globalStiffnessMatrix[3, 2] = globalStiffnessMatrix[2, 3];
            globalStiffnessMatrix[3, 3] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            globalStiffnessMatrix[3, 4] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            globalStiffnessMatrix[3, 5] = 6 * E * I * sinb / (L0 * Lc);

            globalStiffnessMatrix[4, 0] = globalStiffnessMatrix[0, 4];
            globalStiffnessMatrix[4, 1] = globalStiffnessMatrix[1, 4];
            globalStiffnessMatrix[4, 2] = globalStiffnessMatrix[2, 4];
            globalStiffnessMatrix[4, 3] = globalStiffnessMatrix[3, 4];
            globalStiffnessMatrix[4, 4] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[4, 5] = -6 * E * I * cosb / (L0 * Lc);

            globalStiffnessMatrix[5, 0] = globalStiffnessMatrix[0, 5];
            globalStiffnessMatrix[5, 1] = globalStiffnessMatrix[1, 5];
            globalStiffnessMatrix[5, 2] = globalStiffnessMatrix[2, 5];
            globalStiffnessMatrix[5, 3] = globalStiffnessMatrix[3, 5];
            globalStiffnessMatrix[5, 4] = globalStiffnessMatrix[4, 5];
            globalStiffnessMatrix[5, 5] = 4 * E * I / L0;

            return globalStiffnessMatrix;
        }
        #endregion

        public override double[,] CreateLambdaMatrix()
        {
            return null;
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

            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            globalStiffnessMatrix = CreateGlobalStiffnessMatrix();
        }

        public override void CalculateCurrentValues()
        {
            CalculateCurrentNodalCoordinates();
            lengthCurrent = CalculateElementLength(node1XYCurrent, node2XYCurrent);
            sinCurrent = CalculateElementSinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
            cosCurrent = CalculateElementCosinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
            betaAngleCurrent = CalculateElementBetaAngle(node1XYCurrent, node2XYCurrent);
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            globalStiffnessMatrix = CreateGlobalStiffnessMatrix();
        }
    }
}
