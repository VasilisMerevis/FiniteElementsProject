using System;
using System.Collections.Generic;
using System.IO;


namespace FiniteElementsProject
{
	public class InputData
	{
        string[] rawDataLine;
        string[] stringSeparators = new string[] { " " };

        List<string> splitLine;

        public double[] nodesX;
        public double[] nodesY;
        public int[] localnode1;
        public int[] localnode2;
        public string[] elementType;
        public double[] area;
        public double[] elasticity;
        public double[] inertia;
        public int[] boundaryDof;
        public double[] externalForcesVector;
        

        public InputData()
        {
            rawDataLine = System.IO.File.ReadAllLines(@"WriteLines2.txt");
        }

        public void SetNewInputFile(string newFileName)
        {
            rawDataLine = System.IO.File.ReadAllLines(@newFileName);
        }

        private void SplitRawDataLine(int row)
        {
            splitLine = new List<string>(rawDataLine[row].Split(stringSeparators, StringSplitOptions.RemoveEmptyEntries));
        }

        //public void ReadAllData()
        //{
        //    for (int i = 0; i < 9; i++)
        //    {
        //        SplitRawDataLine(i);
        //        splitLine.RemoveAt(0);
        //    }
        //}

        private void ReadNodalCoordinates()
        {
            SplitRawDataLine(0);
            splitLine.RemoveAt(0);
            nodesX = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                nodesX[row] = double.Parse(splitLine[row]);
            }
            SplitRawDataLine(1);
            splitLine.RemoveAt(0);
            nodesY = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                nodesY[row] = Convert.ToDouble(splitLine[row]);
            }
        }

        private void ReadLocalNodes()
        {
            SplitRawDataLine(2);
            splitLine.RemoveAt(0);
            localnode1 = new int[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                localnode1[row] = Convert.ToInt32(splitLine[row]);
            }
            SplitRawDataLine(3);
            splitLine.RemoveAt(0);
            localnode2 = new int[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                localnode2[row] = Convert.ToInt32(splitLine[row]);
            }
        }

        private void ReadElemType()
        {
            SplitRawDataLine(4);
            splitLine.RemoveAt(0);
            elementType = new string[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                elementType[row] = splitLine[row];
            }
        }

        private void ReadElemArea()
        {
            SplitRawDataLine(5);
            splitLine.RemoveAt(0);
            area = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                area[row] = Convert.ToDouble(splitLine[row]);
            }
        }

        private void ReadElemYoungMod()
        {
            SplitRawDataLine(6);
            splitLine.RemoveAt(0);
            elasticity = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                elasticity[row] = Convert.ToDouble(splitLine[row]);
            }
        }

        private void ReadElemInertia()
        {
            SplitRawDataLine(7);
            splitLine.RemoveAt(0);
            inertia = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                inertia[row] = Convert.ToDouble(splitLine[row]);
            }
        }

        private void ReadBoundedDOF()
        {
            SplitRawDataLine(8);
            splitLine.RemoveAt(0);
            boundaryDof = new int[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                boundaryDof[row] = Convert.ToInt32(splitLine[row]);
            }
        }

        private void ReadExternalForces()
        {
            SplitRawDataLine(9);
            splitLine.RemoveAt(0);
            externalForcesVector = new double[splitLine.Count];
            for (int row = 0; row < splitLine.Count; row++)
            {
                externalForcesVector[row] = Convert.ToDouble(splitLine[row]);
            }
        }

        public void ReadAllData()
        {
            ReadNodalCoordinates();
            ReadLocalNodes();
            ReadElemType();
            ReadElemArea();
            ReadElemYoungMod();
            ReadElemInertia();
            ReadBoundedDOF();
            ReadExternalForces();
        }

    }
}

