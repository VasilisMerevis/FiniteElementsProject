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
        public double[] internalLocalForcesVector;
        private double deformationGradient,elasticityTensor, lambda, mi, cauchyStress;

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
            lambda = 0;
            mi = E / 2;
        }

        public override double[,] CreateLambdaMatrix()
        {
            double[,] transformationMatrix = new double[,] { { cosCurrent, sinCurrent, 0, 0 }, { 0, 0, cosCurrent, sinCurrent } };
            return transformationMatrix;
        }

        private double[] TransformToLocalCoordinates(double[] node1XY, double[] node2XY, double cosinus, double sinus)
        {
            double[] globalCoordinatesVector = new double[] { node1XY[0], node1XY[1], node2XY[0], node2XY[1] };
            double[,] transformationMatrix = new double[,] { { cosinus, sinus, 0, 0 }, { 0, 0, cosinus, sinus } };
            double[] localCoordinatesVector = VectorOperations.MatrixVectorProduct(transformationMatrix, globalCoordinatesVector);
            return localCoordinatesVector;
        }

        private double CalculateDeformationGradient()
        {
            double deformationGradient = (currentLocalCoordinatesVector[1] - currentLocalCoordinatesVector[0]) / (initialLocalCoordinatesVector[1] - initialLocalCoordinatesVector[0]);
            return deformationGradient;
        }

        private double CalculateStVenantElasticityTensor()
        {
            double elasticityTensor = lambda + 2 * mi;
            return elasticityTensor;
        }

        private double CalculateNeoHookeElasticityTensor()
        {
            double elasticityTensor = (lambda / deformationGradient) + 2 * (mi - lambda * Math.Log(deformationGradient)) / deformationGradient;
            return elasticityTensor;
        }

        private double CalculateCauchyStressForStVenant()
        {
            double cauchyStress = deformationGradient * ((lambda / 2) * (deformationGradient * deformationGradient - 1) + mi * (deformationGradient * deformationGradient - 1));
            return cauchyStress;
        }

        private double CalculateCauchyStressForNeoHooke()
        {
            double cauchyStress = mi / deformationGradient * ((deformationGradient * deformationGradient - 1) + lambda / deformationGradient * (deformationGradient * deformationGradient - 1));
            return cauchyStress;
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double[,] constitutiveComponentMatrix = new double[,]
            { {-elasticityTensor/Math.Pow(lengthCurrent,2), 0 }, {0, -elasticityTensor/Math.Pow(lengthCurrent,2) } };

            double[,] initialStressComponentmatrix = new double[,]
                {{cauchyStress/Math.Pow(lengthCurrent,2), 0 }, {0, cauchyStress/Math.Pow(lengthCurrent,2) } };

            double[,] fullTangentMatrix = MatrixOperations.ScalarMatrixProduct(A / lengthCurrent, MatrixOperations.MatrixAddition(constitutiveComponentMatrix, initialStressComponentmatrix));
            return fullTangentMatrix;
        }

        public override double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] lambdaTransposeMatrix = MatrixOperations.Transpose(lambdaMatrix);
            double[,] localStiffByLambda = MatrixOperations.MatrixProduct(localStiffnessMatrix, lambdaMatrix);
            double[,] globalStiffnessMatrix4x4 = MatrixOperations.MatrixProduct(lambdaTransposeMatrix, localStiffByLambda);
            globalStiffnessMatrix = new double[,]
                { { globalStiffnessMatrix4x4[0, 0], globalStiffnessMatrix4x4[0,1], 0, globalStiffnessMatrix4x4[0,2], globalStiffnessMatrix4x4[0,3], 0},
                  { globalStiffnessMatrix4x4[1, 0], globalStiffnessMatrix4x4[1,1], 0, globalStiffnessMatrix4x4[1,2], globalStiffnessMatrix4x4[1,3], 0},
                  { 0, 0, 0, 0, 0, 0 },
                  { globalStiffnessMatrix4x4[2, 0], globalStiffnessMatrix4x4[2,1], 0, globalStiffnessMatrix4x4[2,2], globalStiffnessMatrix4x4[2,3], 0},
                  { globalStiffnessMatrix4x4[3, 0], globalStiffnessMatrix4x4[3,1], 0, globalStiffnessMatrix4x4[3,2], globalStiffnessMatrix4x4[3,3], 0},
                  { 0, 0, 0, 0, 0, 0 }
                };

            return globalStiffnessMatrix;
        }

        private double[] CalculateInternalLocalForcesVector()
        {
            double[] internalLocalForcesVector = new double[] {-cauchyStress*A, cauchyStress*A};
            return internalLocalForcesVector;
        }

        private double[] CalculateInternalGlobalForcesVector()
        {
            double[,] transformationMatrix = new double[,] { { cosCurrent, sinCurrent, 0, 0 }, { 0, 0, cosCurrent, sinCurrent } };
            double[] internalGlobalForcesVector4x1 = VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(transformationMatrix), internalLocalForcesVector);
            internalGlobalForcesVector = new double[] { internalGlobalForcesVector4x1[0], internalGlobalForcesVector4x1[1], 0, internalGlobalForcesVector4x1[2], internalGlobalForcesVector4x1[3], 0 };
            return internalGlobalForcesVector;
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
            elasticityTensor = CalculateNeoHookeElasticityTensor();
            cauchyStress = CalculateCauchyStressForNeoHooke();
            //localDisplacementVector = CalculateLocalDisplacementVector();
            internalLocalForcesVector = CalculateInternalLocalForcesVector();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            localStiffnessMatrix = CreateLocalStiffnessMatrix();
            lambdaMatrix = CreateLambdaMatrix();
            globalStiffnessMatrix = CreateGlobalStiffnessMatrix();

        }

        public override void CalculateCurrentValues()
        {
            CalculateCurrentNodalCoordinates();
            lengthCurrent = CalculateElementLength(node1XYCurrent, node2XYCurrent);
            sinCurrent = CalculateElementSinus(node1XYCurrent, node2XYCurrent, lengthCurrent);
            cosCurrent = CalculateElementCosinus(node1XYCurrent, node2XYCurrent, lengthCurrent);

            currentLocalCoordinatesVector = TransformToLocalCoordinates(node1XYCurrent, node2XYCurrent, cosCurrent, sinCurrent);
            //localDisplacementVector = CalculateLocalDisplacementVector();
            deformationGradient = CalculateDeformationGradient();
            elasticityTensor = CalculateNeoHookeElasticityTensor();
            cauchyStress = CalculateCauchyStressForNeoHooke();
            internalLocalForcesVector = CalculateInternalLocalForcesVector();
            internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
            localStiffnessMatrix = CreateLocalStiffnessMatrix();
            lambdaMatrix = CreateLambdaMatrix();
            globalStiffnessMatrix = localStiffnessMatrix;
        }

    }
}
