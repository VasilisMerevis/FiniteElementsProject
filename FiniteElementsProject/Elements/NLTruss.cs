using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class NLTruss : Element1D
    {
        double[] localDisplacementVector;
        public double[] internalLocalForcesVector, internalGlobalForcesVector;

        public NLTruss(double E, double A, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
            this.internalGlobalForcesVector = new double[6];
        }

        private double[] TransformToLocalCoordinates(double[] node1XY, double[] node2XY, double cosinus, double sinus)
        {
            double[] globalCoordinatesVector = new double[] { node1XY[0], node2XY[0], node1XY[1], node2XY[1] };
            double[,] transformationMatrix = new double[,] { { cosinus, sinus, 0, 0 }, { 0, 0, cosinus, sinus } };
            double[] localCoordinatesVector = VectorOperations.MatrixVectorProduct(transformationMatrix, globalCoordinatesVector);
            return localCoordinatesVector;
        }

        public static double Length(double[] vec2d)
        {
            double length = Math.Sqrt(Math.Pow((vec2d[2] - vec2d[0]), 2) + Math.Pow((vec2d[3] - vec2d[1]), 2));
            return length;
        }

        public static double Cosinus(double[] vec2d, double length)
        {
            double cx = (vec2d[2] - vec2d[0]) / length;
            return cx;
        }

        public static double Sinus(double[] vec2d, double length)
        {
            double sx = (vec2d[3] - vec2d[1]) / length;
            return sx;
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

        public static double F11(double[] vec1dUp, double[] vec1dRef)
        {
            double f11 = (vec1dUp[1] - vec1dUp[0]) / (vec1dRef[1] - vec1dRef[0]);
            return f11;
        }

        public static double StVenant(double lamda, double mi)
        {
            double c = lamda + 2 * mi;
            return c;
        }

        public static double NeoHooke(double lamda, double mi, double f11)
        {
            double c = (lamda / f11) + 2 * (mi - lamda * Math.Log(f11)) / f11;
            return c;
        }

        //public override void CalculateInitialValues()
        //{
        //    lengthInitial = CalculateElementLength(node1XYInitial, node2XYInitial);
        //    sinInitial = CalculateElementSinus(node1XYInitial, node2XYInitial, lengthInitial);
        //    cosInitial = CalculateElementCosinus(node1XYInitial, node2XYInitial, lengthInitial);

        //    lengthCurrent = lengthInitial;
        //    sinCurrent = sinInitial;
        //    cosCurrent = cosInitial;

        //    localDisplacementVector = CalculateLocalDisplacementVector();
        //    internalLocalForcesVector = CalculateInternalLocalForcesVector();
        //    internalGlobalForcesVector = CalculateInternalGlobalForcesVector();
        //    localStiffnessMatrix = CreateLocalStiffnessMatrix();
        //    globalStiffnessMatrix = localStiffnessMatrix;

        //}

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
