using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class NLTruss : Element1D
    {
        double[] localDisplacementVector;
        private double[] initialLocalCoordinatesVector, currentLocalCoordinatesVector;
        public double[] internalLocalForcesVector, internalGlobalForcesVector;
        private double deformationGradient;

        public NLTruss(double E, double A, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
            this.internalGlobalForcesVector = new double[6];
            initialLocalCoordinatesVector = new double[2];
            currentLocalCoordinatesVector = new double[2];
        }

        private double[] TransformToLocalCoordinates(double[] node1XY, double[] node2XY, double cosinus, double sinus)
        {
            double[] globalCoordinatesVector = new double[] { node1XY[0], node2XY[0], node1XY[1], node2XY[1] };
            double[,] transformationMatrix = new double[,] { { cosinus, sinus, 0, 0 }, { 0, 0, cosinus, sinus } };
            double[] localCoordinatesVector = VectorOperations.MatrixVectorProduct(transformationMatrix, globalCoordinatesVector);
            return localCoordinatesVector;
        }

        private double CalculateDeformationGradient()
        {
            double deformationGradient = (currentLocalCoordinatesVector[1] - currentLocalCoordinatesVector[0]) / (initialLocalCoordinatesVector[1] - initialLocalCoordinatesVector[0]);
            return deformationGradient;
        }

        private double CalculateStVenantElasticityTensor(double lamda, double mi)
        {
            double elasticityTensor = lamda + 2 * mi;
            return elasticityTensor;
        }

        private double CalculateNeoHookeElasticityTensor(double lamda, double mi, double deformationGradient)
        {
            double elasticityTensor = (lamda / deformationGradient) + 2 * (mi - lamda * Math.Log(deformationGradient)) / deformationGradient;
            return elasticityTensor;
        }

        public static double N1(double ksi)
        {
            double n1 = (1 - ksi) / 2;
            return n1;
        }

        public static double N2(double ksi)
        {
            double n2 = (1 + ksi) / 2;
            return n2;
        }

        public static double Jacobian(double length)
        {
            double j = length / 2;
            return j;
        }

        public static double N1globalDer(double jacobian, double length)
        {
            double n1globalDer = (2 / length) * (-1 / 2);
            return n1globalDer; 
        }

        public static double N2globalDer(double jacobian, double length)
        {
            double n2globalDer = (2 / length) * (1 / 2);
            return n2globalDer;
        }

        

        

        public override void CalculateInitialValues()
        {
            lengthInitial = CalculateElementLength(node1XYInitial, node2XYInitial);
            sinInitial = CalculateElementSinus(node1XYInitial, node2XYInitial, lengthInitial);
            cosInitial = CalculateElementCosinus(node1XYInitial, node2XYInitial, lengthInitial);

            lengthCurrent = lengthInitial;
            sinCurrent = sinInitial;
            cosCurrent = cosInitial;

            initialLocalCoordinatesVector = TransformToLocalCoordinates(node1XYInitial, node2XYInitial, cosInitial, sinInitial);
            currentLocalCoordinatesVector = initialLocalCoordinatesVector;

            deformationGradient = CalculateDeformationGradient();

            //localDisplacementVector = CalculateLocalDisplacementVector();
            //internalLocalForcesVector = CalculateInternalLocalForcesVector();
            //internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            //localStiffnessMatrix = CreateLocalStiffnessMatrix();
            //globalStiffnessMatrix = localStiffnessMatrix;

        }

        //public override void CalculateCurrentValues()
        //{
        //    CalculateCurrentNodalCoordinates();
        //    lengthCurrent = CalculateElementLength(node1XYCurrent, node2XYCurrent);
        //    sinCurrent = CalculateElementSinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
        //    cosCurrent = CalculateElementCosinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
        //    localDisplacementVector = CalculateLocalDisplacementVector();
        //    internalLocalForcesVector = CalculateInternalLocalForcesVector();
        //    internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
        //    localStiffnessMatrix = CreateLocalStiffnessMatrix();
        //    globalStiffnessMatrix = localStiffnessMatrix;
        //}

    }
}
