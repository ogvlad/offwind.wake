using System;
using System.IO;

namespace WakeCode
{
    public class DataReader
    {
        //----------------------------------------------------
        //************************************************
        //  SUBROUTINE READ THE DATA !
        //------------------------------------------------
        public void READ_DATA(SolverData solverData, GeneralData generalData)
        {
            using (var fileStream = File.Open("initial_data.inp", FileMode.OpenOrCreate, FileAccess.Read))
            using (var streamReader = new StreamReader(fileStream))
            {
                generalData.IMAX = ReadInt(streamReader); // The number of grid points in x direction
                generalData.JMAX = ReadInt(streamReader); // The number of the grid points in Y direction

                generalData.x = new double[generalData.IMAX];
                generalData.y = new double[generalData.JMAX];
                generalData.vell_i = new double[generalData.IMAX,generalData.JMAX];

                solverData.Dturb = ReadDouble(streamReader); // THE DIAMETER OF THE TURBIN
                solverData.H = ReadDouble(streamReader); //  THE HEIGHT OF THE TURBINE
                solverData.Ct = ReadDouble(streamReader); // TURBINE THRUST COEFFICIENT
                solverData.Kwake = ReadDouble(streamReader); // wake expand scalar
                solverData.Uhub = ReadDouble(streamReader);
                    //m/s - VELOCITY AT THE HUB, WITHOUT THE INFLUENCE OF THE WIND TURBIN
                generalData.N_TURB = ReadInt(streamReader); //THE NUMBER OF THE TURBINE

                generalData.x_turb = new double[generalData.N_TURB];
                generalData.y_turb = new double[generalData.N_TURB];
                generalData.R_TURB = new double[generalData.N_TURB];
                generalData.WPOWER = new double[generalData.N_TURB];

                generalData.xc_turb = new System.Int32[generalData.N_TURB];
                generalData.yc_turb = new System.Int32[generalData.N_TURB];

                solverData.Rho = ReadDouble(streamReader); // THE DENSITY OF THE AIR 
                solverData.dist = ReadDouble(streamReader); // the distance behind the turbine where the power is computed
                generalData.ang = ReadDouble(streamReader);
                    // rotational angle of the axis: vellocity has the same direction as Ox
                ReadEmpty(streamReader);
                ReadEmpty(streamReader);
                for (var i = 0; i <= generalData.N_TURB - 1; i++)
                {
                    var t = ReadXY(streamReader); // pozition of the turbine
                    generalData.x_turb[i] = t.Item1;
                    generalData.y_turb[i] = t.Item2;
                }
                ReadEmpty(streamReader);
            }
        }

        private static void ReadEmpty(TextReader textReader)
        {
            textReader.ReadLine();
        }

        private static int ReadInt(TextReader textReader)
        {
            int intValue;
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 1) || !int.TryParse(lineParts[0], out intValue))
            {
                throw new FormatException();
            }
            return intValue;
        }

        private static double ReadDouble(TextReader textReader)
        {
            double doubleValue;
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 1) || !double.TryParse(lineParts[0], out doubleValue))
            {
                throw new FormatException();
            }
            return doubleValue;
        }

        private static Tuple<double, double> ReadXY(TextReader textReader)
        {
            double doubleValue1;
            double doubleValue2;
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 2) || !double.TryParse(lineParts[0], out doubleValue1) || !double.TryParse(lineParts[1], out doubleValue2))
            {
                throw new FormatException();
            }
            return new Tuple<double, double>(doubleValue1, doubleValue2);
        }
    }
}