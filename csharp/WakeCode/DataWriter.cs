﻿using System.Linq;

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

            using (System.IO.FileStream fileStream = System.IO.File.Open("FLOW.xyz", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (var streamWriter = new System.IO.StreamWriter(fileStream))
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
            }

            using (System.IO.FileStream fileStream = System.IO.File.Open("FLOW.q", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (var streamWriter = new System.IO.StreamWriter(fileStream))
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
                            WRITE(streamWriter, solverData.Rho * generalData.vell_i[i, j]);
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
        }

        /// <summary>
        /// SUBROUTINE  _DATA Power
        /// </summary>
        /// <param name="generalData"></param>
        public void WRITE_DATA_power(GeneralData generalData)
        {
            using (System.IO.FileStream fileStream = System.IO.File.Open("Power_Output.dat", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (System.IO.StreamWriter streamWriter = new System.IO.StreamWriter(fileStream))
                {
                    WRITE(streamWriter, "   Turbine Number(m)   ", "Turbine Location-X(m)   ", "Turbine Location-Y(m)    ", "POWER(W)");
                    int i;
                    for (i = 1; i <= generalData.N_TURB; i++)
                    {
                        WRITE(streamWriter, i, generalData.x_turb[i - 1], generalData.y_turb[i - 1], generalData.WPOWER[i - 1]);
                    }
                }
            }
        }

        private void WRITE(System.IO.TextWriter textWriter, params object[] values)
        {
            string line = values.Aggregate("", (string partialLine, object value) => { return partialLine + " " + value.ToString(); }, (string partialLine) => { return (partialLine.Length >= 1 ? partialLine.Substring(1) : partialLine); });
            textWriter.WriteLine(line);
        }
    }
}
