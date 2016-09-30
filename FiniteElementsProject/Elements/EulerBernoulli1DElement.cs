using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    class EulerBernoulli1DElement : Element1D
    {
        private double I;

        public EulerBernoulli1DElement(double E, double A, double I, double[] nodesX, double[] nodesY)
            :base(E, A, nodesX, nodesY)
        {
            this.I = I;
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double length = VectorOperations.CalculateVectorLengthFromEndPoints (nodesX [0], nodesX [1], nodesY [0], nodesY [1]);
            localStiffnessMatrix = new[,]
                { { E * A / length, 0, 0, -E * A / length, 0, 0 },
                { 0, 12 * E * I / Math.Pow(length, 3), 6 * E * I / Math.Pow(length, 2), 0, -12 * E * I / Math.Pow(length, 3), 6 * E * I / Math.Pow( length, 2) },
                { 0, 6 * E * I / Math.Pow(length, 2), 4 * E * I / length, 0, -6 * E * I / Math.Pow(length, 2), 2 * E * I / length},
                { -E * A / length, 0, 0, E* A/ length, 0, 0 },
                { 0, -12 * E * I / Math.Pow(length, 3), -6 * E * I / Math.Pow(length, 2), 0, 12 * E * I / Math.Pow(length, 3), -6 * E * I / Math.Pow(length, 2) },
                { 0, 6 * E * I / Math.Pow( length, 2), 2 * E * I / length, 0, -6 * E * I / Math.Pow(length, 2), 4 * E * I / length}};

            return localStiffnessMatrix;
        }

        public override double[,] CreateLambdaMatrix()
        {
            double length = Math.Sqrt (Math.Pow ((nodesY [1] - nodesY [0]), 2) + Math.Pow ((nodesX [1] - nodesX [0]), 2));
            double cosine = (nodesX [1] - nodesX [0]) / length;
            double sine = (nodesY [1] - nodesY [0]) / length;
            lambdaMatrix = new[,]
            { { cosine, sine, 0, 0, 0, 0 },
            { -sine, cosine, 0, 0, 0, 0 },
            { 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, cosine, sine, 0 },
            { 0, 0, 0, -sine, cosine, 0 },
            { 0, 0, 0, 0, 0, 1 }};

            return lambdaMatrix;
        }
    }
}
