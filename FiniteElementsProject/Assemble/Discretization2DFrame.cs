using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    public class Discretization2DFrame
    {
		private string[] elementType;
		protected double[] E, A, I, nodesX, nodesY;
		protected int[] localNode1, localNode2;
		protected Element1D[] beamElementsList;
		private double[,] totalStiffnessMatrix, totalMassMatrix;
        public double[] internalForcesTotalVector;


        public double[,] TotalStiffnessMatrix {
			get{ return totalStiffnessMatrix;}
		}

        public double[,] TotalMassMatrix
        {
            get { return totalMassMatrix; }
        }

        public Discretization2DFrame(InputData inputData) 
        {
			this.elementType =inputData.elementType;
			this.E = inputData.elasticity;
			this.A = inputData.area;
			this.I = inputData.inertia;
			this.nodesX = inputData.nodesX;
			this.nodesY = inputData.nodesY;
			this.localNode1 = inputData.localnode1;
			this.localNode2 = inputData.localnode2;
            beamElementsList = new Element1D[localNode1.Length];
            totalStiffnessMatrix = new double [3*nodesX.Length,3*nodesX.Length];
            totalMassMatrix = new double[3 * nodesX.Length, 3 * nodesX.Length];
            this.internalForcesTotalVector = new double[3 * nodesX.Length];
        }

        public Element1D[] GetStiffnessMatrices()
        {
            for (int elem = 0; elem < localNode1.Length; elem++)
            {
				double[] elementNodesX = {nodesX[localNode1[elem]-1] , nodesX[localNode2[elem]-1]};
				double[] elementNodesY = {nodesY[localNode1[elem]-1] , nodesY[localNode2[elem]-1]};

				switch(elementType[elem])
				{
                    case "Beam":
                        beamElementsList[elem] = new EulerBernoulli1DElement(E[elem], A[elem], I[elem], elementNodesX, elementNodesY);
                        break;
                    case "Bar":
                        beamElementsList[elem] = new Bar1DElement(E[elem], A[elem], elementNodesX, elementNodesY);
                        break;
                    case "NLBeam":
                        beamElementsList[elem] = new NLEulerBernoulli1DElement(E[elem], A[elem], I[elem], elementNodesX, elementNodesY);
                        break;
                    case "NLTruss":
                        beamElementsList[elem] = new NLTruss(E[elem], A[elem], elementNodesX, elementNodesY);
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

        private int[] ElementDOFs(int elementNumber)
        {
            int[] dof = new int[6];
            dof[0] = localNode1[elementNumber] * 3 - 2;
            dof[1] = localNode1[elementNumber] * 3 - 1;
            dof[2] = localNode1[elementNumber] * 3;

            dof[3] = localNode2[elementNumber] * 3 - 2;
            dof[4] = localNode2[elementNumber] * 3 - 1;
            dof[5] = localNode2[elementNumber] * 3;
            return dof;
        }

        public void UpdateValues(double[] totalDisplacementVector)
        {
            int totalElements = localNode1.Length;
            for (int element = 0; element < totalElements; element++)
            {
                int node1 = localNode1[element];
                int node2 = localNode2[element];
                double[] node1GlobalDisplacementVector = new double[] { totalDisplacementVector[node1 * 3 - 2-1], totalDisplacementVector[node1 * 3 - 1-1], totalDisplacementVector[node1 * 3-1] };
                double[] node2GlobalDisplacementVector = new double[] { totalDisplacementVector[node2 * 3 - 2-1], totalDisplacementVector[node2 * 3 - 1-1], totalDisplacementVector[node2 * 3-1] };
                beamElementsList[element].SetGlobalDisplacementVector(node1GlobalDisplacementVector, node2GlobalDisplacementVector);
                beamElementsList[element].CalculateCurrentValues();
            }
        }

        public double[,] CreateTotalStiffnessMatrix()
        {
            for (int element = 0; element < localNode1.Length; element++)
            {
                int[] dof = ElementDOFs(element);

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        totalStiffnessMatrix[dof[i] - 1, dof[j] - 1] = totalStiffnessMatrix[dof[i] - 1, dof[j] - 1] + beamElementsList[element].globalStiffnessMatrix[i, j];
                    }
                }
            }
            return totalStiffnessMatrix;
        }

        public double[,] CreateTotalMassMatrix()
        {
            for (int element = 0; element < localNode1.Length; element++)
            {
                int[] dof = ElementDOFs(element);

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        totalMassMatrix[dof[i] - 1, dof[j] - 1] = totalMassMatrix[dof[i] - 1, dof[j] - 1] + beamElementsList[element].massMatrix[i, j];
                    }
                }
            }
            return totalMassMatrix;
        }

        public double[] CreateTotalInternalForcesVector()
        {
            int totalNodes = nodesX.Length;
            int totalElements = localNode1.Length;

            for (int element = 0; element < totalElements; element++)
            {
                int[] dof = { localNode1[element] * 3 - 2, localNode1[element] * 3 - 1, localNode1[element]*3, localNode2[element] * 3 - 2, localNode2[element] * 3 - 1, localNode2[element]*3 };
                for (int i = 0; i < 6; i++)
                {
                    internalForcesTotalVector[dof[i] - 1] = internalForcesTotalVector[dof[i] - 1] + beamElementsList[element].internalGlobalForcesVector[i];
                }
            }
            return internalForcesTotalVector;
        }
    }
}
