using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    public class Assembly : IAssembly
    {
        private string[] elementType;
        protected double[] E, A, I, nodesX, nodesY;
        protected int[] localNode1, localNode2;
        protected Element1D[] beamElementsList;
        private double[,] totalStiffnessMatrix, totalMassMatrix;
        public double[] internalForcesTotalVector;
        private int totalDOF;
        private int[] boundedDOFsVector;



        public double[,] TotalStiffnessMatrix
        {
            get { return totalStiffnessMatrix; }
        }

        public double[,] TotalMassMatrix
        {
            get { return totalMassMatrix; }
        }

        public int[] BoundedDOFsVector
        {
            get
            {
                return boundedDOFsVector;
            }

            set
            {
                boundedDOFsVector = value;
            }
        }

        public Assembly(InputData inputData)
        {
            this.elementType = inputData.elementType;
            this.E = inputData.elasticity;
            this.A = inputData.area;
            this.I = inputData.inertia;
            this.nodesX = inputData.nodesX;
            this.nodesY = inputData.nodesY;
            this.localNode1 = inputData.localnode1;
            this.localNode2 = inputData.localnode2;
            beamElementsList = new Element1D[localNode1.Length];
            totalStiffnessMatrix = new double[3 * nodesX.Length, 3 * nodesX.Length];
            totalMassMatrix = new double[3 * nodesX.Length, 3 * nodesX.Length];
            this.internalForcesTotalVector = new double[3 * nodesX.Length];
        }

        public void InitializeMatrices()
        {
            this.totalStiffnessMatrix = new double[totalDOF, totalDOF];
            this.totalMassMatrix = new double[totalDOF, totalDOF];
            this.internalForcesTotalVector = new double[totalDOF];
        }

        public Element1D[] GetStiffnessMatrices()
        {
            for (int elem = 0; elem < localNode1.Length; elem++)
            {
                double[] elementNodesX = { nodesX[localNode1[elem] - 1], nodesX[localNode2[elem] - 1] };
                double[] elementNodesY = { nodesY[localNode1[elem] - 1], nodesY[localNode2[elem] - 1] };

                switch (elementType[elem])
                {
                    case "Beam":
                        beamElementsList[elem] = new EulerBernoulli1DElement(E[elem], A[elem], I[elem], elementNodesX, elementNodesY);
                        totalDOF = 3 * nodesX.Length;
                        break;
                    case "Bar":
                        beamElementsList[elem] = new Bar1DElement(E[elem], A[elem], elementNodesX, elementNodesY);
                        totalDOF = 2 * nodesX.Length;
                        break;
                    case "Bar2D":
                        beamElementsList[elem] = new Bar2D(E[elem], A[elem], elementNodesX, elementNodesY);
                        totalDOF = 2 * nodesX.Length;
                        break;
                    case "NLBeam":
                        beamElementsList[elem] = new NLEulerBernoulli1DElement(E[elem], A[elem], I[elem], elementNodesX, elementNodesY);
                        totalDOF = 3 * nodesX.Length;
                        break;
                    case "NLTruss":
                        beamElementsList[elem] = new NLTruss(E[elem], A[elem], elementNodesX, elementNodesY);
                        totalDOF = 3 * nodesX.Length;
                        break;
                    case "Contact":
                        beamElementsList[elem] = new ContactNTN2D(E[elem], A[elem], elementNodesX, elementNodesY);
                        totalDOF = 3 * nodesX.Length;
                        break;

                }

                beamElementsList[elem].CalculateInitialValues();
            }
            return beamElementsList;
        }

        public void GetMassMatrices()
        {
            for (int elem = 0; elem < localNode1.Length; elem++)
            {
                beamElementsList[elem].CreateMassMatrix();
            }
        }

        public void UpdateValues(double[] totalDisplacementVector)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(totalDisplacementVector, boundedDOFsVector);
            int totalElements = localNode1.Length;
            for (int element = 0; element < totalElements; element++)
            {
                int node1 = localNode1[element];
                int node2 = localNode2[element];
                double[] node1GlobalDisplacementVector = new double[] { fullTotalDisplacementVector[node1 * 3 - 2 - 1], fullTotalDisplacementVector[node1 * 3 - 1 - 1], fullTotalDisplacementVector[node1 * 3 - 1] };
                double[] node2GlobalDisplacementVector = new double[] { fullTotalDisplacementVector[node2 * 3 - 2 - 1], fullTotalDisplacementVector[node2 * 3 - 1 - 1], fullTotalDisplacementVector[node2 * 3 - 1] };
                beamElementsList[element].SetGlobalDisplacementVector(node1GlobalDisplacementVector, node2GlobalDisplacementVector);
                beamElementsList[element].CalculateCurrentValues();
            }
        }

        public double[,] CreateTotalStiffnessMatrix()
        {
            Array.Clear(TotalStiffnessMatrix, 0, TotalStiffnessMatrix.Length);
            for (int element = 0; element < localNode1.Length; element++)
            {
                List<int> dof = beamElementsList[element].ElementDOFs(localNode1, localNode2, element);
                for (int i = 0; i < dof.Count; i++)
                {
                    for (int j = 0; j < dof.Count; j++)
                    {
                        totalStiffnessMatrix[dof[i] - 1, dof[j] - 1] = totalStiffnessMatrix[dof[i] - 1, dof[j] - 1] + beamElementsList[element].globalStiffnessMatrix[i, j];
                    }
                }
            }
            double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalStiffnessMatrix, boundedDOFsVector);
            return reducedStiffnessMatrix;
        }

        public double[,] CreateTotalMassMatrix()
        {
            Array.Clear(TotalMassMatrix, 0, TotalMassMatrix.Length);
            for (int element = 0; element < localNode1.Length; element++)
            {
                List<int> dof = beamElementsList[element].ElementDOFs(localNode1, localNode2, element);
                for (int i = 0; i < dof.Count; i++)
                {
                    for (int j = 0; j < dof.Count; j++)
                    {
                        totalMassMatrix[dof[i] - 1, dof[j] - 1] = totalMassMatrix[dof[i] - 1, dof[j] - 1] + beamElementsList[element].massMatrix[i, j];
                    }
                }
            }
            double[,] reducedMassMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalMassMatrix, boundedDOFsVector);
            return reducedMassMatrix;
        }

        public double[] CreateTotalInternalForcesVector()
        {
            Array.Clear(internalForcesTotalVector, 0, internalForcesTotalVector.Length);
            int totalNodes = nodesX.Length;
            int totalElements = localNode1.Length;

            for (int element = 0; element < totalElements; element++)
            {
                int[] dof = { localNode1[element] * 3 - 2, localNode1[element] * 3 - 1, localNode1[element] * 3, localNode2[element] * 3 - 2, localNode2[element] * 3 - 1, localNode2[element] * 3 };
                for (int i = 0; i < 6; i++)
                {
                    internalForcesTotalVector[dof[i] - 1] = internalForcesTotalVector[dof[i] - 1] + beamElementsList[element].internalGlobalForcesVector[i];
                }
            }
            double[] reducedInternalForcesVector = BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundedDOFsVector);
            return reducedInternalForcesVector;
        }
    }
}
