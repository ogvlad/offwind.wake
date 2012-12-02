using System;

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
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int i, j;

            using (System.IO.FileStream fileStream = System.IO.File.Open("initial_data.inp", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Read))
            {
                using (System.IO.StreamReader streamReader = new System.IO.StreamReader(fileStream))
                {
                    READ(streamReader, ref generalData.IMAX);               // The number of grid points in x direction
                    READ(streamReader, ref generalData.JMAX);               // The number of the grid points in Y direction

                    generalData.x = new double[generalData.IMAX];
                    generalData.y = new double[generalData.JMAX];
                    generalData.vell_i = new double[generalData.IMAX, generalData.JMAX];

                    READ(streamReader, ref solverData.Dturb);               // THE DIAMETER OF THE TURBIN
                    READ(streamReader, ref solverData.H);                   //  THE HEIGHT OF THE TURBINE
                    READ(streamReader, ref solverData.Ct);                  // TURBINE THRUST COEFFICIENT
                    READ(streamReader, ref solverData.Kwake);               // wake expand scalar
                    READ(streamReader, ref solverData.Uhub);                //m/s - VELOCITY AT THE HUB, WITHOUT THE INFLUENCE OF THE WIND TURBIN
                    READ(streamReader, ref generalData.N_TURB);             //THE NUMBER OF THE TURBINE

                    generalData.x_turb = new double[generalData.N_TURB];
                    generalData.y_turb = new double[generalData.N_TURB];
                    generalData.R_TURB = new double[generalData.N_TURB];
                    generalData.WPOWER = new double[generalData.N_TURB];

                    generalData.xc_turb = new System.Int32[generalData.N_TURB];
                    generalData.yc_turb = new System.Int32[generalData.N_TURB];

                    READ(streamReader, ref solverData.Rho);                         // THE DENSITY OF THE AIR 
                    READ(streamReader, ref solverData.dist);                        // the distance behind the turbine where the power is computed
                    READ(streamReader, ref generalData.ang);                        // rotational angle of the axis: vellocity has the same direction as Ox
                    READ(streamReader);
                    READ(streamReader);
                    for (i = 0; i <= generalData.N_TURB - 1; i++)
                    {
                        READ(streamReader, ref generalData.x_turb[i], ref generalData.y_turb[i]);   // pozition of the turbine
                    }
                    READ(streamReader);
                }
            }
        }  // END SUBROUTINE READ DATA


        private static void READ(System.IO.TextReader textReader)
        {
            textReader.ReadLine();
        }

        private static void READ(System.IO.TextReader textReader, ref int intValue)
        {
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 1) || !int.TryParse(lineParts[0], out intValue))
            {
                throw new FormatException();
            }
        }

        private static void READ(System.IO.TextReader textReader, ref double doubleValue)
        {
            double result;
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 1) || !double.TryParse(lineParts[0], out result))
            {
                throw new FormatException();
            }
            doubleValue = result;
        }

        private static void READ(System.IO.TextReader textReader, ref double doubleValue1, ref double doubleValue2)
        {
            string line = textReader.ReadLine();
            string[] lineParts = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (!(lineParts.Length >= 2) || !double.TryParse(lineParts[0], out doubleValue1) || !double.TryParse(lineParts[1], out doubleValue2))
            {
                throw new FormatException();
            }
        }
    }
}