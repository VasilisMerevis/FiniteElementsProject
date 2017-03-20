using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class ContactNTN2D : Element1D
    {
        private double PenaltyFactor { get; set; }
        private double[] normalUnitVector;
        private double[,] oldTangentMatrix;
        private double penetration;
        public double[,] TangentMatrix { get; set; }
        public ContactNTN2D(double E, double A, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            lambdaMatrix = new double[4, 4];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[4, 4];
            oldTangentMatrix = new double[4, 4];
        }

        private double[] CalculateNormalUnitVector()
        {
            double[] normalVector = new double[] { nodesX[1] - nodesX[0], nodesY[1] - nodesY[0] };
            double normalVectorLength = VectorOperations.VectorNorm2(normalVector);
            double[] normalUnitVec = new double[] { normalVector[0] / normalVectorLength, normalVector[1] / normalVectorLength };
            return normalUnitVec;
        }

        private double[,] CalculatePositionMatrix()
        {
            double[,] aMatrix = new double[,]
                {
                    { -1,0,1,0},
                    {0,-1,0,1 }
                };
            return aMatrix;
        }

        private double CalculateKsi3()
        {
            double[,] A = CalculatePositionMatrix();
            double[,] AT = MatrixOperations.Transpose(A);
            double [] n = CalculateNormalUnitVector();
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = new double[] { node1XYCurrent[0], node1XYCurrent[1], node2XYCurrent[0], node2XYCurrent[1] };
            double ksi3 = VectorOperations.VectorDotProduct(xupd, AT_n);
            return ksi3;
        }

        public override double[,] CreateGlobalStiffnessMatrix()
        {
            double[] n = CalculateNormalUnitVector();
            double[,] A = CalculatePositionMatrix();
            double[,] AT = MatrixOperations.Transpose(A);
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
            double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
            double[,] AT_nxn_A = MatrixOperations.MatrixProduct(AT, nxn_A);
            double[,] e_AT_nxn_A = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
            double[,] globalStiffnessMatrix = MatrixOperations.MatrixAddition(oldTangentMatrix, e_AT_nxn_A);
            return globalStiffnessMatrix;
        }

        public double[] CalculateInternalGlobalForcesVector()
        {
            double[,] A = CalculatePositionMatrix();
            double[,] AT = MatrixOperations.Transpose(A);
            double[] n = CalculateNormalUnitVector();
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double ksi = CalculateKsi3();
            double[] ksi_AT_n = VectorOperations.VectorScalarProductNew(AT_n, ksi);
            double[] e_ksi_AT_n = VectorOperations.VectorScalarProductNew(ksi_AT_n, PenaltyFactor);
            return e_ksi_AT_n;
        }

        public override void CalculateInitialValues()
        {
            penetration = 0.0;
            globalStiffnessMatrix = CreateGlobalStiffnessMatrix();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
        }

        public override void CalculateCurrentValues()
        {
            CalculateCurrentNodalCoordinates();
            penetration = CalculateKsi3();
            if (penetration <= 0)
            {
                globalStiffnessMatrix = MatrixOperations.MatrixAddition(globalStiffnessMatrix, CreateGlobalStiffnessMatrix());
                internalGlobalForcesVector = VectorOperations.VectorVectorSubtraction(internalGlobalForcesVector, CalculateInternalGlobalForcesVector());
            }
        }
    }
}
