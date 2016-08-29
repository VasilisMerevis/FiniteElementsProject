using System;

namespace FiniteElementsProject
{
	public interface IElements1D
	{
		double[,] globalStiffnessMatrix { get; set; }
		double[,] CreateLocalStiffnessMatrix ();
		double[,] CreateLambdaMatrix ();
		double[,] CreateGlobalStiffnessMatrix ();
	}
}

