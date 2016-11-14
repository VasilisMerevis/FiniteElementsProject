using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class Bar2D : Element1D
    {

        public double density = 8000;

        public Bar2D(double E, double A, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            lambdaMatrix = new double[4, 4];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[4, 4];
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double length = VectorOperations.CalculateVectorLengthFromEndPoints(nodesX[0], nodesX[1], nodesY[0], nodesY[1]);
            localStiffnessMatrix = new[,]
            { { E * A / length, 0, 0, -E * A / length, 0, 0 },
                { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0},
                { -E * A / length, 0, 0, E* A/ length, 0, 0 },
                { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0}};

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

        public double[,] CreateGlobalStiffnessmatrix()
        {
            double L = Math.Sqrt(Math.Pow((nodesY[1] - nodesY[0]), 2) + Math.Pow((nodesX[1] - nodesX[0]), 2));
            double c = (nodesX[1] - nodesX[0]) / L;
            double s = (nodesY[1] - nodesY[0]) / L;
            double cs = c * s;
            double c2 = c * c;
            double s2 = s * s;
            double[,] globalStiffnessMatrix = new double[,]
                {
                    {A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
                    {A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
                    {-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
                    {-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
                };
            return globalStiffnessMatrix;
        }

        public override double[,] CreateMassMatrix()
        {
            double length = Math.Sqrt(Math.Pow((nodesY[1] - nodesY[0]), 2) + Math.Pow((nodesX[1] - nodesX[0]), 2));
            double elementMass = density * A * length;
            massMatrix = new double[,]
            {
                {elementMass/2, 0, 0, 0},
                {0, elementMass/2, 0, 0},
                {0, 0, elementMass/2, 0},
                {0, 0, 0, elementMass/2},
            };

            return massMatrix;
        }

        public override void CalculateInitialValues()
        {
            globalStiffnessMatrix = this.CreateGlobalStiffnessmatrix();
        }
    }
}
