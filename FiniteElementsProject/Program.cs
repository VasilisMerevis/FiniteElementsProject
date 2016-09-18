using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
	class Program
	{
        static void Main(string[] args)
        {
            //Declaration of basic frame info
            InputData data = new InputData();
            data.ReadAllData();

            if (data.elementType.Contains("NLBeam"))// | data.elementType.Contains("NLTruss"))
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();
                //Exercise1Frame.CreateTotalStiffnessMatrix();
                NLSolver solution = new NLSolver(Exercise1Frame, 1000, data);
                solution.SolveWithMethod("Load Controlled Newton-Raphson");
                VectorOperations.PrintVector(solution.solutionVector);
            }
            else if (data.elementType.Contains("NLTruss"))
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();

                //Exercise1Frame.CreateTotalStiffnessMatrix();
                NLSolver solution = new NLSolver(Exercise1Frame, 10000, data);
                solution.SolveWithMethod("Newton-Raphson");
                VectorOperations.PrintVector(solution.solutionVector);
            }
            else
            {
                //Creation of local, global and total stiffness matrices
                Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
                Exercise1Frame.GetStiffnessMatrices();
                Exercise1Frame.CreateTotalStiffnessMatrix();

                //Creation reduced matrix depended on boundary conditions
                double[,] reducedTotalStiff = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalStiffnessMatrix, data.boundaryDof);

                //Solution using Cholesky factorization with forward and backward substitution
                DirectSolver solution = new DirectSolver(reducedTotalStiff, data.externalForcesVector);
                solution.SolveWithMethod("Gauss");

                VectorOperations.PrintVector(solution.GetSolutionVector);

                Console.WriteLine();
                IterativeSolver solution2 = new IterativeSolver(reducedTotalStiff, data.externalForcesVector, 1000);
                solution2.SolveWithMethod("PCG");

                VectorOperations.PrintVector(solution2.GetSolutionVector);
            }
            
        }
    }
}
