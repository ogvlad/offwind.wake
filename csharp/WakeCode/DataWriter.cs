using System.IO;
using System.Linq;

namespace WakeCode
{
    public class DataWriter
    {
        /// <summary>
        /// SUBROUTINE  _DATA
        /// </summary>
        /// <param name="solverData"></param>
        /// <param name="generalData"></param>
        public void WRITE_DATA(SolverData solverData, GeneralData generalData)
        {
            int i;
            int j;

            using (FileStream fileStream = File.Open("FLOW.xyz", FileMode.OpenOrCreate, FileAccess.Write))
            using (var streamWriter = new StreamWriter(fileStream))
            {
                WRITE(streamWriter, generalData.IMAX, generalData.JMAX);
                for (j = 1; j <= generalData.JMAX; j++)
                {
                    for (i = 0; i <= generalData.IMAX - 1; i++)
                    {
                        WRITE(streamWriter, generalData.x[i]);
                    }
                }
                WRITE(streamWriter);
                for (j = 0; j <= generalData.JMAX - 1; j++)
                {
                    for (i = 1; i <= generalData.IMAX; i++)
                    {
                        WRITE(streamWriter, generalData.y[j]);
                    }
                }
            }

            using (FileStream fileStream = File.Open("FLOW.q", FileMode.OpenOrCreate, FileAccess.Write))
            using (var streamWriter = new StreamWriter(fileStream))
            {
                WRITE(streamWriter, generalData.IMAX, generalData.JMAX);
                WRITE(streamWriter, "0.1   ", "  10  ", "  10000  ", "  0.1 ");
                for (j = 1; j <= generalData.JMAX; j++)
                {
                    for (i = 1; i <= generalData.IMAX; i++)
                    {
                        WRITE(streamWriter, solverData.Rho);
                    }
                }
                for (j = 0; j <= generalData.JMAX - 1; j++)
                {
                    for (i = 0; i <= generalData.IMAX - 1; i++)
                    {
                        WRITE(streamWriter, solverData.Rho*generalData.vell_i[i, j]);
                    }
                }
                for (j = 1; j <= generalData.JMAX; j++)
                {
                    for (i = 1; i <= generalData.IMAX; i++)
                    {
                        WRITE(streamWriter, 0);
                    }
                }
                for (j = 1; j <= generalData.JMAX; j++)
                {
                    for (i = 1; i <= generalData.IMAX; i++)
                    {
                        WRITE(streamWriter, 0);
                    }
                }
            }
        }

        /// <summary>
        /// SUBROUTINE  _DATA Power
        /// </summary>
        /// <param name="generalData"></param>
        public void WRITE_DATA_power(GeneralData generalData)
        {
            using (var fileStream = File.Open("Power_Output.dat", FileMode.OpenOrCreate, FileAccess.Write))
            using (var streamWriter = new StreamWriter(fileStream))
            {
                WRITE(streamWriter, "   Turbine Number(m)   ", "Turbine Location-X(m)   ",
                      "Turbine Location-Y(m)    ", "POWER(W)");
                int i;
                for (i = 1; i <= generalData.N_TURB; i++)
                {
                    WRITE(streamWriter, i, generalData.x_turb[i - 1], generalData.y_turb[i - 1],
                          generalData.WPOWER[i - 1]);
                }
            }
        }

        private void WRITE(TextWriter textWriter, params object[] values)
        {
            string line = values.Aggregate("",
                                           (partialLine, value) => partialLine + " " + value,
                                           partialLine => (partialLine.Length >= 1
                                                               ? partialLine.Substring(1)
                                                               : partialLine));
            textWriter.WriteLine(line);
        }
    }
}