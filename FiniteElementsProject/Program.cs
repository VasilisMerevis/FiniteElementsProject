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
            ////Declaration of basic frame info
            //InputData data = new InputData ();
            //         data.ReadAllData();

            ////Creation of local, global and total stiffness matrices
            //Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
            //         Exercise1Frame.GetStiffnessMatrices();
            //         Exercise1Frame.CreateTotalStiffnessMatrix();

            //         //Creation reduced matrix depended on boundary conditions
            //         double[,] reducedTotalStiff = BoundaryConditionsImposition.ReducedTotalStiff(Exercise1Frame.TotalStiffnessMatrix, data.boundaryDof);

            //         //Solution using Cholesky factorization with forward and backward substitution
            //         DirectSolver solution = new DirectSolver(reducedTotalStiff, data.externalForcesVector);
            //         solution.SolveWithMethod("Gauss");

            //         VectorOperations.PrintVector(solution.GetSolutionVector);

            //         Console.WriteLine();
            //         IterativeSolver solution2 = new IterativeSolver(reducedTotalStiff, data.externalForcesVector, 1000);
            //         solution2.SolveWithMethod("PCG");

            //         VectorOperations.PrintVector(solution2.GetSolutionVector);

            //         //IterativeSolver solution3 = new IterativeSolver(reducedTotalStiff, data.externalForcesVector, 1000);
            //         //solution3.SolveWithMethod("Jacobi");

            //         //VectorOperations.PrintVector(solution3.GetSolutionVector);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Non Linear
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Declaration of basic frame info
            InputData data = new InputData();
            data.ReadAllData();

            //Creation of local, global and total stiffness matrices
            Discretization2DFrame Exercise1Frame = new Discretization2DFrame(data);
            Exercise1Frame.GetStiffnessMatrices();
            //Exercise1Frame.CreateTotalStiffnessMatrix();

            NLSolver solution = new NLSolver(Exercise1Frame, 1000, data);


            solution.SolveWithMethod("Newton-Raphson");
            VectorOperations.PrintVector(solution.solutionVector);
        }
    }
}
