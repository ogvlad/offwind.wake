using System;
using System.Linq;

namespace WakeCode
{
    public class WakeCalc
    {
        public const double pi = 3.1415926535897;

        public void Run(GeneralData generalData, SolverData solverData)
        {
            double ppp;

            //************************************************************************
            //ROTATE THE DOMAIN, AND THE X,Y COORDINATE OF THE TURBINE so that the wind to be in x direction
            //------------------------------------------------------------------
            ROTATE_coord(generalData);
            if (solverData.Ct > 1)
            {
                Console.WriteLine(" The value of the Ct should be less 1, hence Ct=0.3)");
                solverData.Ct = 0.3;
            }

            solverData.Cp = 0.5 * (1 + Math.Sqrt(1 - solverData.Ct)) * solverData.Ct;

            ORDER(generalData);

            ppp = 5.0;

            DOMAIN_PT(ref generalData.x, ref generalData.IMAX, ref generalData.dx, ref solverData.Dturb, ref generalData.x_turb, ref generalData.N_TURB, ref generalData.xmax, ref generalData.xmin, ref ppp);

            ppp = 2.0;

            DOMAIN_PT(ref generalData.y, ref generalData.JMAX, ref generalData.dy, ref solverData.Dturb, ref generalData.y_turb, ref generalData.N_TURB, ref generalData.ymax, ref generalData.ymin, ref ppp);
            Turb_centr_coord(ref generalData.N_TURB, ref generalData.IMAX, ref generalData.x, ref generalData.x_turb, ref generalData.xc_turb);
            Turb_centr_coord(ref generalData.N_TURB, ref generalData.JMAX, ref generalData.y, ref generalData.y_turb, ref generalData.yc_turb);
            COMPUTE_VELL(solverData, generalData);

            WRITE_DATA(solverData, generalData);

            COMPUTE_WPower(solverData, generalData);
            WRITE_DATA_power(generalData);
        }

        /// <summary>
        /// rotate the coordinate of the turbines
        /// </summary>
        /// <param name="generalData"></param>
        public void ROTATE_coord(GeneralData generalData)
        {
            int i;
            var XX_TURB = new double[generalData.N_TURB];
            var YY_TURB = new double[generalData.N_TURB];
            double ang1;
            ang1 = generalData.ang * pi / 180;
            for (i = 0; i <= generalData.N_TURB - 1; i++)
            {
                XX_TURB[i] = generalData.x_turb[i] * Math.Cos(ang1) - generalData.y_turb[i] * Math.Sin(ang1);
                YY_TURB[i] = generalData.x_turb[i] * Math.Sin(ang1) + generalData.y_turb[i] * Math.Cos(ang1);
            }

            for (i = 0; i <= generalData.N_TURB - 1; i++)
            {
                generalData.x_turb[i] = XX_TURB[i];
                generalData.y_turb[i] = YY_TURB[i];
            }
        } // 


        /// <summary>
        /// COMPUTE THE GRID POINTS
        /// </summary>
        /// <param name="XX"></param>
        /// <param name="IIMAX"></param>
        /// <param name="DDX"></param>
        /// <param name="DDtur"></param>
        /// <param name="XX_TURB"></param>
        /// <param name="NN_TURB"></param>
        /// <param name="XXMAX"></param>
        /// <param name="XXMIN"></param>
        /// <param name="pppoint"></param>
        private void DOMAIN_PT(ref double[] XX, ref System.Int32 IIMAX, ref double DDX, ref double DDtur, ref double[] XX_TURB, ref System.Int32 NN_TURB, ref double XXMAX, ref double XXMIN, ref double pppoint)
        {
            int i;

            XXMAX = XX_TURB[0];
            XXMIN = XX_TURB[0];
            for (i = 1; i <= NN_TURB - 1; i++)
            {
                if (XX_TURB[i] > XXMAX)
                {
                    XXMAX = XX_TURB[i];
                }
                if (XX_TURB[i] < XXMIN)
                {
                    XXMIN = XX_TURB[i];
                }
            }

            XXMAX = XXMAX + DDtur * pppoint;
            XXMIN = XXMIN - 2 * DDtur;

            XX[0] = XXMIN;
            DDX = (XXMAX - XXMIN) / (IIMAX - 1);
            for (i = 1; i <= IIMAX - 1; i++)
            {
                XX[i] = XX[i - 1] + DDX;
            }
            pppoint = 0.0;
        }

        /// <summary>
        /// ROTATE COORDINATE THAT THE wind to be in X ? DIRECTION
        /// the subroutine determine the coordinates of the center of the turbine
        /// </summary>
        /// <param name="nn"></param>
        /// <param name="iimax"></param>
        /// <param name="xx"></param>
        /// <param name="xx_turb"></param>
        /// <param name="xxc_turb"></param>
        private void Turb_centr_coord(ref System.Int32 nn, ref System.Int32 iimax, ref double[] xx, ref double[] xx_turb, ref System.Int32[] xxc_turb)
        {
            System.Int32 i;
            System.Int32 ii;

            for (i = 0; i <= nn - 1; i++)
            {
                for (ii = 0; ii <= iimax - 2; ii++)
                {
                    if (xx[ii] <= xx_turb[i] && xx_turb[i] < xx[ii + 1])
                    {
                        xxc_turb[i] = ii + 1;
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// ORDER of THE TURBINE in function of x coordinate
        /// </summary>
        /// <param name="generalData"></param>
        private void ORDER(GeneralData generalData)
        {
            int i;
            int j;
            int k;
            double aa;
            double bb;

            for (i = 1; i <= generalData.N_TURB - 1; i++)
            {
                for (j = 0; j <= generalData.N_TURB - i - 1; j++)
                {
                    if (generalData.x_turb[j] > generalData.x_turb[j + 1])
                    {
                        aa = generalData.x_turb[j];
                        bb = generalData.y_turb[j];
                        generalData.x_turb[j] = generalData.x_turb[j + 1];
                        generalData.y_turb[j] = generalData.y_turb[j + 1];
                        generalData.x_turb[j + 1] = aa;
                        generalData.y_turb[j + 1] = bb;
                    }
                }
            }

            //for (j=1; j <= N_turb; j++) {
            for (i = 0; i <= generalData.N_TURB - 1; i++)
            {
                for (k = i + 1; k <= generalData.N_TURB - 1; k++)
                {
                    if (generalData.x_turb[i] == generalData.x_turb[k])
                    {
                        if (generalData.y_turb[i] > generalData.y_turb[k])
                        {
                            aa = generalData.y_turb[i];
                            generalData.y_turb[i] = generalData.y_turb[k];
                            generalData.y_turb[k] = aa;
                        }
                    }
                }
            }
            //}

        }

        private void WRITE(System.IO.TextWriter textWriter, params object[] values)
        {
            string line = values.Aggregate("", (string partialLine, object value) => { return partialLine + " " + value.ToString(); }, (string partialLine) => { return (partialLine.Length >= 1 ? partialLine.Substring(1) : partialLine); });
            textWriter.WriteLine(line);
        }

        /// <summary>
        /// SUBROUTINE  _DATA
        /// </summary>
        /// <param name="solverData"></param>
        /// <param name="generalData"></param>
        private void WRITE_DATA(SolverData solverData, GeneralData generalData)
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

        private int INT(double doubleValue)
        {
            return (int)doubleValue;
        }

        /// <summary>
        /// SUBROUTINE  _SHADOW AREA
        /// </summary>
        /// <param name="solverData"></param>
        /// <param name="generalData"></param>
        private void COMPUTE_VELL(SolverData solverData, GeneralData generalData)
        {
            int I;
            int J;
            int K;
            int Ni;
            int nj;
            int jj_max;
            int jj_min;
            int ii;
            int nk;
            double[] SHADOW = new double[generalData.N_TURB];
            double DIJ;
            double RR_I;
            double ALPHA_I;
            double ALPHA_K;
            double LIJ;
            double PP;
            double SS;
            double ss0;
            double RR_k;
            double vv;

            double r0;
            double x_dist;
            double rr_max = 0;
            double rrt;
            double area = 0;

            for (I = 0; I <= generalData.IMAX - 1; I++)
            {
                for (J = 0; J <= generalData.JMAX - 1; J++)
                {
                    generalData.vell_i[I, J] = solverData.Uhub;
                }
            }

            r0 = 0.5 * solverData.Dturb; // all the tubine have the same diameter

            nk = 2 * (INT(solverData.Dturb / generalData.dy));
            for (K = 1; K <= generalData.N_TURB; K++)
            {
                J = 0;
                SS = 0.0;
                ss0 = (pi * r0 * r0);

                for (I = 1; I <= K - 1; I = I + 1) // calculate the influence of the turbine i over the turbine k
                {
                    RR_I = r0 + solverData.Kwake * (generalData.x[generalData.xc_turb[K - 1] - 1] - generalData.x[generalData.xc_turb[I - 1] - 1]);
                    DIJ = Math.Abs(generalData.y_turb[I - 1] - generalData.y_turb[K - 1]);
                    if (RR_I >= (r0 + DIJ) || DIJ <= generalData.dy)
                    {
                        SS = SS + ((r0 * r0) / (RR_I * RR_I));
                    }
                    else
                    {
                        if ((DIJ) < (RR_I + r0) && (DIJ) > generalData.dy)
                        {
                            J = J + 1;
                            ALPHA_I = (RR_I * RR_I) + (DIJ * DIJ) - (r0 * r0);
                            ALPHA_I = ALPHA_I / (2 * RR_I * DIJ);
                            ALPHA_I = Math.Acos(ALPHA_I);
                            ALPHA_K = (r0 * r0) + (DIJ * DIJ) - (RR_I * RR_I);
                            ALPHA_K = ALPHA_K / (2 * r0 * DIJ);
                            ALPHA_K = Math.Acos(ALPHA_K);
                            AAREA(ref RR_I, ref r0, ref DIJ, ref area);

                            SHADOW[J - 1] = (ALPHA_I * (Math.Pow(RR_I, 2)) + ALPHA_K * (Math.Pow(r0, 2)));
                            SHADOW[J - 1] = SHADOW[J - 1] - 2 * area;
                            SS = SS + ((SHADOW[J - 1]) / ss0) * ((r0 * r0) / (RR_I * RR_I));
                        }
                        else
                        {
                            SS = SS;
                        }
                    }
                }

                for (ii = generalData.xc_turb[K - 1]; ii <= generalData.IMAX; ii++)
                {
                    rrt = r0 + solverData.Kwake * (generalData.x[ii - 1] - generalData.x[generalData.xc_turb[K - 1] - 1]);
                    rr_max = Math.Max(rrt, rr_max);
                    nj = (INT(rrt / generalData.dy));
                    jj_min = Math.Max(1, generalData.yc_turb[K - 1] - nj);
                    jj_max = Math.Min(generalData.JMAX, generalData.yc_turb[K - 1] + nj);

                    for (J = jj_min; J <= jj_max; J++)
                    {
                        if (((-generalData.vell_i[ii - 1, J - 1] + solverData.Uhub) > 0) && (ii > generalData.xc_turb[K - 1] + nk))
                        {
                            vv = generalData.vell_i[ii - 1, J - 1];
                            generalData.vell_i[ii - 1, J - 1] = solverData.Uhub + solverData.Uhub * (Math.Sqrt(1 - solverData.Ct) - 1) * ((r0 * r0) / (rrt * rrt));
                            generalData.vell_i[ii - 1, J - 1] = generalData.vell_i[ii - 1, J - 1] * (1 - (1 - Math.Sqrt(1 - solverData.Ct)) * SS);
                            //vell_i(ii,j)=(vell_i(ii,j)+0.15*vv)/1.15;
                            generalData.vell_i[ii - 1, J - 1] = Math.Min(vv, generalData.vell_i[ii - 1, J - 1]);
                        }
                        else
                        {
                            generalData.vell_i[ii - 1, J - 1] = solverData.Uhub + solverData.Uhub * (Math.Sqrt(1 - solverData.Ct) - 1) * (r0 / rrt) * (r0 / rrt);
                            generalData.vell_i[ii - 1, J - 1] = generalData.vell_i[ii - 1, J - 1] * (1 - (1 - Math.Sqrt(1 - solverData.Ct)) * SS);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// SUBROUTINE  compute the power at the distance dist behind the WT
        /// </summary>
        /// <param name="solverData"></param>
        /// <param name="generalData"></param>
        private void COMPUTE_WPower(SolverData solverData, GeneralData generalData)
        {
            int I;
            int J;
            int K;
            int Ni;
            int nj;
            int jj_max;
            int jj_min;
            int ii;
            int JJ;
            int nd;
            int nk;
            int mm;
            double[] SHADOW = new double[generalData.N_TURB];
            double DIJ;
            double RR_I;
            double ALPHA_I;
            double ALPHA_K;
            double LIJ;
            double PP;
            double SS;
            double ss0;
            double RR_k;
            double vv;
            double SPOWER;
            double vv1;
            double vv2;
            double[] v_power = new double[generalData.N_TURB];

            double r0;
            double x_dist;
            double rr_max;
            double rrt;
            double area = 0;

            r0 = 0.5 * solverData.Dturb; // all the tubine have the same diameter

            ss0 = (pi * r0 * r0);
            I = (int)Math.Truncate(solverData.dist / generalData.dx);
            nd = Math.Max(1, I);

            for (K = 1; K <= generalData.N_TURB; K++)
            {
                J = 0;
                SS = 0.0;
                nk = Math.Max(1, generalData.xc_turb[K - 1] - nd);
                vv1 = generalData.vell_i[nk - 1, generalData.yc_turb[K - 1] - 1];
                generalData.WPOWER[K - 1] = 0.0;
                vv2 = 0.0;
                for (I = K - 1; I >= 1; I = I - 1) // calculate the influence of the turbine i over the turbine k
                {
                    RR_I = r0 + solverData.Kwake * (generalData.x[nk - 1] - generalData.x[generalData.xc_turb[I - 1] - 1]);
                    DIJ = Math.Abs(generalData.y_turb[I - 1] - generalData.y_turb[K - 1]);

                    if (((DIJ) < (RR_I + r0)) && (RR_I <= (r0 + DIJ)))
                    {
                        J = J + 1;
                        ALPHA_I = (RR_I * RR_I) + (DIJ * DIJ) - (r0 * r0);
                        ALPHA_I = ALPHA_I / (2 * RR_I * DIJ);
                        ALPHA_I = Math.Acos(ALPHA_I);
                        ALPHA_K = (r0 * r0) + (DIJ * DIJ) - (RR_I * RR_I);
                        ALPHA_K = ALPHA_K / (2 * r0 * DIJ);
                        ALPHA_K = Math.Acos(ALPHA_K);
                        AAREA(ref RR_I, ref r0, ref DIJ, ref area);

                        SHADOW[J - 1] = (ALPHA_I * (Math.Pow(RR_I, 2)) + ALPHA_K * (Math.Pow(r0, 2)));
                        SHADOW[J - 1] = SHADOW[J - 1] - 2 * area;

                        SS = SS + SHADOW[J - 1];
                        if (SS < ss0)
                        {
                            if (generalData.y_turb[K - 1] > generalData.y_turb[I - 1])
                            {
                                mm = (INT(RR_I / generalData.dy));
                                jj_max = Math.Min(generalData.JMAX, generalData.yc_turb[I - 1] + mm + 1);
                                jj_min = Math.Max(1, generalData.yc_turb[I - 1] + mm - 2);
                                vv1 = generalData.vell_i[nk - 1, jj_max - 1];
                                v_power[J - 1] = generalData.vell_i[nk - 1, jj_min - 1];
                            }
                            else
                            {
                                mm = (INT(RR_I / generalData.dy));
                                jj_max = Math.Min(generalData.JMAX, generalData.yc_turb[I - 1] + mm + 1);
                                jj_min = Math.Max(1, generalData.yc_turb[I - 1] + mm - 2);
                                vv1 = generalData.vell_i[nk - 1, jj_min - 1];
                                v_power[J - 1] = generalData.vell_i[nk - 1, jj_max - 1];
                            }
                        }
                        else
                        {
                            J = J - 1;
                        }
                    }
                }

                if (J > 0)
                {
                    for (I = 1; I <= J; I++)
                    {
                        vv2 = v_power[J - 1] * SHADOW[J - 1] + vv2;
                    }
                }
                vv2 = (vv2 + vv1 * (ss0 - SS)) / ss0;

                generalData.WPOWER[K - 1] = 0.5 * solverData.Rho * (Math.Pow(vv2, 3)) * ss0 * solverData.Cp;
            }
        }

        /// <summary>
        /// SUBROUTINE  _DATA Power
        /// </summary>
        /// <param name="generalData"></param>
        private void WRITE_DATA_power(GeneralData generalData)
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

        /// <summary>
        /// FUNCTION : COMPUTE AREA
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <param name="area"></param>
        private void AAREA(ref double X, ref double Y, ref double Z, ref double area)
        {
            double PP = (X + Y + Z) * 0.5;
            area = Math.Sqrt(PP * (PP - X) * (PP - Y) * (PP - Z));
            return;
        }
    }
}
