using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace WakeCode
{
    internal class GeneralData
    {
        public static System.Int32 N_TURB;
        public static System.Int32 IMAX;
        public static System.Int32 JMAX;
        public static double dx;
        public static double dy;
        public static double pi;
        public static double xmax;
        public static double ymax;
        public static double ymin;
        public static double xmin;
        public static double ang;

        public static double[] x;
        public static double[] y;
        public static double[,] vell_i;
        public static double[] x_turb;     // location of the turbine
        public static double[] y_turb;     // location of the turbine
        public static double[] R_TURB;     // location of the turbine
        public static double[] WPOWER;     // location of the turbine
        public static Int32[] xc_turb;
        public static Int32[] yc_turb;
    }

    internal class SolverData
    {
        public static double Ct;
        public static double Dturb;
        public static double Kwake;
        public static double H;
        public static double Uhub;
        public static double Dwake;
        public static double Rho, Cp, dist;

        public static float[,] V = new float[1000, 1000];
        public static float[,] Darea = new float[1000, 1000];
        public static float[,] Darea_D = new float[1000, 1000];
    }

    static class Program
    {
        private const double pi = 3.1415926535897;
        //private const double pi = Math.PI;

        private static void Main(string[] args)
        {
            //****** declaration of the variable *******************************
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            System.Int32 i, j, k;
            double ppp;
            float W_coef;
            char[] file1 = new char[80];
            //*************************************************

            //************************************************************************
            //ROTATE THE DOMAIN, AND THE X,Y COORDINATE OF THE TURBINE so that the wind to be in x direction
            //------------------------------------------------------------------
            READ_DATA();
            ROTATE_coord();
            if (SolverData.Ct > 1)
            {
                Console.WriteLine(" The value of the Ct should be less 1, hence Ct=0.3)");
                SolverData.Ct = 0.3;
            }

            SolverData.Cp = 0.5 * (1 + Math.Sqrt(1 - SolverData.Ct)) * SolverData.Ct;

            ORDER();

            ppp = 5.0;

            DOMAIN_PT(ref GeneralData.x, ref GeneralData.IMAX, ref GeneralData.dx, ref SolverData.Dturb, ref GeneralData.x_turb, ref GeneralData.N_TURB, ref GeneralData.xmax, ref GeneralData.xmin, ref ppp);

            ppp = 2.0;

            DOMAIN_PT(ref GeneralData.y, ref GeneralData.JMAX, ref GeneralData.dy, ref SolverData.Dturb, ref GeneralData.y_turb, ref GeneralData.N_TURB, ref GeneralData.ymax, ref GeneralData.ymin, ref ppp);
            Turb_centr_coord(ref GeneralData.N_TURB, ref GeneralData.IMAX, ref GeneralData.x, ref GeneralData.x_turb, ref GeneralData.xc_turb);
            Turb_centr_coord(ref GeneralData.N_TURB, ref GeneralData.JMAX, ref GeneralData.y, ref GeneralData.y_turb, ref GeneralData.yc_turb);
            COMPUTE_VELL();

            WRITE_DATA();

            COMPUTE_WPower();
            WRITE_DATA_power();
        }

        //*********************************************************
        //   rotate the coordinate of the turbine
        //--------------------------------------------
        private static void ROTATE_coord()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData;

            int i;
            double[] XX_TURB = new double[GeneralData.N_TURB], YY_TURB = new double[GeneralData.N_TURB];
            double ang1;
            ang1 = GeneralData.ang * pi / 180;
            for (i = 0; i <= GeneralData.N_TURB - 1; i++)
            {
                XX_TURB[i] = GeneralData.x_turb[i] * Math.Cos(ang1) - GeneralData.y_turb[i] * Math.Sin(ang1);
                YY_TURB[i] = GeneralData.x_turb[i] * Math.Sin(ang1) + GeneralData.y_turb[i] * Math.Cos(ang1);
            }

            for (i = 0; i <= GeneralData.N_TURB - 1; i++)
            {
                GeneralData.x_turb[i] = XX_TURB[i];
                GeneralData.y_turb[i] = YY_TURB[i];
            }
        } // 

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

        //----------------------------------------------------
        //************************************************
        //  SUBROUTINE READ THE DATA !
        //------------------------------------------------
        private static void READ_DATA()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int i, j;

            using (System.IO.FileStream fileStream = System.IO.File.Open("initial_data.inp", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Read))
            {
                using (System.IO.StreamReader streamReader = new System.IO.StreamReader(fileStream))
                {
                    READ(streamReader, ref GeneralData.IMAX);               // The number of grid points in x direction
                    READ(streamReader, ref GeneralData.JMAX);               // The number of the grid points in Y direction

                    GeneralData.x = new double[GeneralData.IMAX];
                    GeneralData.y = new double[GeneralData.JMAX];
                    GeneralData.vell_i = new double[GeneralData.IMAX, GeneralData.JMAX];

                    READ(streamReader, ref SolverData.Dturb);               // THE DIAMETER OF THE TURBIN
                    READ(streamReader, ref SolverData.H);                   //  THE HEIGHT OF THE TURBINE
                    READ(streamReader, ref SolverData.Ct);                  // TURBINE THRUST COEFFICIENT
                    READ(streamReader, ref SolverData.Kwake);               // wake expand scalar
                    READ(streamReader, ref SolverData.Uhub);                //m/s - VELOCITY AT THE HUB, WITHOUT THE INFLUENCE OF THE WIND TURBIN
                    READ(streamReader, ref GeneralData.N_TURB);             //THE NUMBER OF THE TURBINE

                    GeneralData.x_turb = new double[GeneralData.N_TURB];
                    GeneralData.y_turb = new double[GeneralData.N_TURB];
                    GeneralData.R_TURB = new double[GeneralData.N_TURB];
                    GeneralData.WPOWER = new double[GeneralData.N_TURB];

                    GeneralData.xc_turb = new System.Int32[GeneralData.N_TURB];
                    GeneralData.yc_turb = new System.Int32[GeneralData.N_TURB];

                    READ(streamReader, ref SolverData.Rho);                         // THE DENSITY OF THE AIR 
                    READ(streamReader, ref SolverData.dist);                        // the distance behind the turbine where the power is computed
                    READ(streamReader, ref GeneralData.ang);                        // rotational angle of the axis: vellocity has the same direction as Ox
                    READ(streamReader);
                    READ(streamReader);
                    for (i = 0; i <= GeneralData.N_TURB - 1; i++)
                    {
                        READ(streamReader, ref GeneralData.x_turb[i], ref GeneralData.y_turb[i]);   // pozition of the turbine
                    }
                    READ(streamReader);
                }
            }
        }  // END SUBROUTINE READ DATA


        //--------------------------------------------------------------------
        //   COMPUTE THE GRID POINTS 
        //************************************************************************
        private static void DOMAIN_PT(ref double[] XX, ref System.Int32 IIMAX, ref double DDX, ref double DDtur, ref double[] XX_TURB, ref System.Int32 NN_TURB, ref double XXMAX, ref double XXMIN, ref double pppoint)
        {
            int i;
            //REAL(kind=8) ::XX(IIMAX)
            //REAL (kind=8) ::XX_TURB(NN_TURB) 

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
        } // subroutine that compute the grid points
        //-----------------------------------------------------------


        //************************************************************?
        //   ROTATE COORDINATE THAT THE wind to be in X ? DIRECTION
        //------------------------------------------------------------- 

        //***********************************************************************
        // the subroutine determine the coordinates of the center of the turbine
        //*****************************************************
        private static void Turb_centr_coord(ref System.Int32 nn, ref System.Int32 iimax, ref double[] xx, ref double[] xx_turb, ref System.Int32[] xxc_turb)
        {
            System.Int32 i, ii;
            //INTEGER (kind=4):: xxc_turb(nn)
            //REAL (kind=8):: xx_turb(nn), xx(iimax)

            for (i = 0; i <= nn - 1; i++)
            {
                for (ii = 0; ii <= iimax - 2; ii++)
                {
                    if (xx[ii] <= xx_turb[i] && xx_turb[i] < xx[ii + 1])
                    {
                        xxc_turb[i] = ii + 1;
                        //EXIT
                        break;
                    }
                }
            }
        } // subroutine deternines the center coordinates of the turbines
        //--------------------------------------------------------------

        //*****************************************************************
        // ORDER of THE TURBINE in function of x coordinate
        //_******************************************************************??? 
        private static void ORDER()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int i, j, k;
            double aa, bb;

            for (i = 1; i <= GeneralData.N_TURB - 1; i++)
            {
                for (j = 0; j <= GeneralData.N_TURB - i - 1; j++)
                {
                    if (GeneralData.x_turb[j] > GeneralData.x_turb[j + 1])
                    {
                        aa = GeneralData.x_turb[j];
                        bb = GeneralData.y_turb[j];
                        GeneralData.x_turb[j] = GeneralData.x_turb[j + 1];
                        GeneralData.y_turb[j] = GeneralData.y_turb[j + 1];
                        GeneralData.x_turb[j + 1] = aa;
                        GeneralData.y_turb[j + 1] = bb;
                    }
                }
            }

            //for (j=1; j <= N_turb; j++) {
            for (i = 0; i <= GeneralData.N_TURB - 1; i++)
            {
                for (k = i + 1; k <= GeneralData.N_TURB - 1; k++)
                {
                    if (GeneralData.x_turb[i] == GeneralData.x_turb[k])
                    {
                        if (GeneralData.y_turb[i] > GeneralData.y_turb[k])
                        {
                            aa = GeneralData.y_turb[i];
                            GeneralData.y_turb[i] = GeneralData.y_turb[k];
                            GeneralData.y_turb[k] = aa;
                        }
                    }
                }
            }
            //}

        } // SUBROUTINE puts in order the turbine
        //--------------------------------------------------------------------------

        private static void WRITE(System.IO.TextWriter textWriter, params object[] values)
        {
            string line = values.Aggregate("", (string partialLine, object value) => { return partialLine + " " + value.ToString(); }, (string partialLine) => { return (partialLine.Length >= 1 ? partialLine.Substring(1) : partialLine); });
            textWriter.WriteLine(line);
        }

        //************************************************************************
        // *                         SUBROUTINE  _DATA                        *
        //*********************************************************************** 
        private static void WRITE_DATA()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int i, j;
            //	character*80 file1
            //	file1='datab.1'

            using (System.IO.FileStream fileStream = System.IO.File.Open("FLOW.xyz", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (System.IO.StreamWriter streamWriter = new System.IO.StreamWriter(fileStream))
                {
                    WRITE(streamWriter, GeneralData.IMAX, GeneralData.JMAX);
                    for (j = 1; j <= GeneralData.JMAX; j++)
                    {
                        for (i = 0; i <= GeneralData.IMAX - 1; i++)
                        {
                            WRITE(streamWriter, GeneralData.x[i]);
                        }
                    }
                    WRITE(streamWriter);
                    for (j = 0; j <= GeneralData.JMAX - 1; j++)
                    {
                        for (i = 1; i <= GeneralData.IMAX; i++)
                        {
                            WRITE(streamWriter, GeneralData.y[j]);
                        }
                    }
                }
            }

            using (System.IO.FileStream fileStream = System.IO.File.Open("FLOW.q", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (System.IO.StreamWriter streamWriter = new System.IO.StreamWriter(fileStream))
                {
                    WRITE(streamWriter, GeneralData.IMAX, GeneralData.JMAX);
                    WRITE(streamWriter, "0.1   ", "  10  ", "  10000  ", "  0.1 ");
                    for (j = 1; j <= GeneralData.JMAX; j++)
                    {
                        for (i = 1; i <= GeneralData.IMAX; i++)
                        {
                            WRITE(streamWriter, SolverData.Rho);
                        }
                    }
                    for (j = 0; j <= GeneralData.JMAX - 1; j++)
                    {
                        for (i = 0; i <= GeneralData.IMAX - 1; i++)
                        {
                            WRITE(streamWriter, SolverData.Rho * GeneralData.vell_i[i, j]);
                        }
                    }
                    for (j = 1; j <= GeneralData.JMAX; j++)
                    {
                        for (i = 1; i <= GeneralData.IMAX; i++)
                        {
                            WRITE(streamWriter, 0);
                        }
                    }
                    for (j = 1; j <= GeneralData.JMAX; j++)
                    {
                        for (i = 1; i <= GeneralData.IMAX; i++)
                        {
                            WRITE(streamWriter, 0);
                        }
                    }
                }
            }
        }
        //************************************************************************

        private static int INT(double doubleValue)
        {
            return (int)doubleValue;
        }

        //************************************************************************
        // *                         SUBROUTINE  _SHADOW AREA 
        //*********************************************************************** 
        private static void COMPUTE_VELL()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int I, J, K, Ni, nj, jj_max, jj_min, ii, nk;
            double[] SHADOW = new double[GeneralData.N_TURB];
            double DIJ, RR_I, ALPHA_I, ALPHA_K, LIJ;
            double PP, SS, ss0, RR_k, vv;

            double r0, x_dist, rr_max = 0, rrt, area = 0;

            for (I = 0; I <= GeneralData.IMAX - 1; I++)
            {
                for (J = 0; J <= GeneralData.JMAX - 1; J++)
                {
                    GeneralData.vell_i[I, J] = SolverData.Uhub;
                }
            }

            r0 = 0.5 * SolverData.Dturb; // all the tubine have the same diameter

            nk = 2 * (INT(SolverData.Dturb / GeneralData.dy));
            for (K = 1; K <= GeneralData.N_TURB; K++)
            {
                J = 0;
                SS = 0.0;
                ss0 = (pi * r0 * r0);

                for (I = 1; I <= K - 1; I = I + 1) // calculate the influence of the turbine i over the turbine k
                {
                    RR_I = r0 + SolverData.Kwake * (GeneralData.x[GeneralData.xc_turb[K - 1] - 1] - GeneralData.x[GeneralData.xc_turb[I - 1] - 1]);
                    DIJ = Math.Abs(GeneralData.y_turb[I - 1] - GeneralData.y_turb[K - 1]);
                    if (RR_I >= (r0 + DIJ) || DIJ <= GeneralData.dy)
                    {
                        SS = SS + ((r0 * r0) / (RR_I * RR_I));
                    }
                    else
                    {
                        if ((DIJ) < (RR_I + r0) && (DIJ) > GeneralData.dy)
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

                for (ii = GeneralData.xc_turb[K - 1]; ii <= GeneralData.IMAX; ii++)
                {
                    rrt = r0 + SolverData.Kwake * (GeneralData.x[ii - 1] - GeneralData.x[GeneralData.xc_turb[K - 1] - 1]);
                    rr_max = Math.Max(rrt, rr_max);
                    nj = (INT(rrt / GeneralData.dy));
                    jj_min = Math.Max(1, GeneralData.yc_turb[K - 1] - nj);
                    jj_max = Math.Min(GeneralData.JMAX, GeneralData.yc_turb[K - 1] + nj);

                    for (J = jj_min; J <= jj_max; J++)
                    {
                        if (((-GeneralData.vell_i[ii - 1, J - 1] + SolverData.Uhub) > 0) && (ii > GeneralData.xc_turb[K - 1] + nk))
                        {
                            vv = GeneralData.vell_i[ii - 1, J - 1];
                            GeneralData.vell_i[ii - 1, J - 1] = SolverData.Uhub + SolverData.Uhub * (Math.Sqrt(1 - SolverData.Ct) - 1) * ((r0 * r0) / (rrt * rrt));
                            GeneralData.vell_i[ii - 1, J - 1] = GeneralData.vell_i[ii - 1, J - 1] * (1 - (1 - Math.Sqrt(1 - SolverData.Ct)) * SS);
                            //vell_i(ii,j)=(vell_i(ii,j)+0.15*vv)/1.15;
                            GeneralData.vell_i[ii - 1, J - 1] = Math.Min(vv, GeneralData.vell_i[ii - 1, J - 1]);
                        }
                        else
                        {
                            GeneralData.vell_i[ii - 1, J - 1] = SolverData.Uhub + SolverData.Uhub * (Math.Sqrt(1 - SolverData.Ct) - 1) * (r0 / rrt) * (r0 / rrt);
                            GeneralData.vell_i[ii - 1, J - 1] = GeneralData.vell_i[ii - 1, J - 1] * (1 - (1 - Math.Sqrt(1 - SolverData.Ct)) * SS);
                        }
                    }
                }
            }
        }   // subroutine that compute the velocity in front of the wind turbine


        //******************************************************************************

        //
        //************************************************************************
        // *       SUBROUTINE  compute the power at the distance dist behind the WT
        //*********************************************************************** 
        private static void COMPUTE_WPower()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int I, J, K, Ni, nj, jj_max, jj_min, ii, JJ, nd, nk, mm;
            double[] SHADOW = new double[GeneralData.N_TURB];
            double DIJ, RR_I, ALPHA_I, ALPHA_K, LIJ;
            double PP, SS, ss0, RR_k, vv, SPOWER, vv1, vv2;
            double[] v_power = new double[GeneralData.N_TURB];

            double r0, x_dist, rr_max, rrt, area = 0;

            r0 = 0.5 * SolverData.Dturb; // all the tubine have the same diameter

            ss0 = (pi * r0 * r0);
            I = (int)Math.Truncate(SolverData.dist / GeneralData.dx);
            nd = Math.Max(1, I);

            for (K = 1; K <= GeneralData.N_TURB; K++)
            {
                J = 0;
                SS = 0.0;
                nk = Math.Max(1, GeneralData.xc_turb[K - 1] - nd);
                vv1 = GeneralData.vell_i[nk - 1, GeneralData.yc_turb[K - 1] - 1];
                GeneralData.WPOWER[K - 1] = 0.0;
                vv2 = 0.0;
                for (I = K - 1; I >= 1; I = I - 1) // calculate the influence of the turbine i over the turbine k
                {
                    RR_I = r0 + SolverData.Kwake * (GeneralData.x[nk - 1] - GeneralData.x[GeneralData.xc_turb[I - 1] - 1]);
                    DIJ = Math.Abs(GeneralData.y_turb[I - 1] - GeneralData.y_turb[K - 1]);

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
                            if (GeneralData.y_turb[K - 1] > GeneralData.y_turb[I - 1])
                            {
                                mm = (INT(RR_I / GeneralData.dy));
                                jj_max = Math.Min(GeneralData.JMAX, GeneralData.yc_turb[I - 1] + mm + 1);
                                jj_min = Math.Max(1, GeneralData.yc_turb[I - 1] + mm - 2);
                                vv1 = GeneralData.vell_i[nk - 1, jj_max - 1];
                                v_power[J - 1] = GeneralData.vell_i[nk - 1, jj_min - 1];
                            }
                            else
                            {
                                mm = (INT(RR_I / GeneralData.dy));
                                jj_max = Math.Min(GeneralData.JMAX, GeneralData.yc_turb[I - 1] + mm + 1);
                                jj_min = Math.Max(1, GeneralData.yc_turb[I - 1] + mm - 2);
                                vv1 = GeneralData.vell_i[nk - 1, jj_min - 1];
                                v_power[J - 1] = GeneralData.vell_i[nk - 1, jj_max - 1];
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

                GeneralData.WPOWER[K - 1] = 0.5 * SolverData.Rho * (Math.Pow(vv2, 3)) * ss0 * SolverData.Cp;
            }
        }   // subroutine that compute the velocity in front of the wind turbine

        //**********************************************************************************

        //************************************************************************
        // *                         SUBROUTINE  _DATA Power                       *
        //*********************************************************************** 
        private static void WRITE_DATA_power()
        {
            //GeneralData GeneralData = new GeneralData();
            //SolverData SolverData = new SolverData();
            int i;
            //char[] file1 = new char[80];
            //file1 = "power_data.1";

            using (System.IO.FileStream fileStream = System.IO.File.Open("Power_Output.dat", System.IO.FileMode.OpenOrCreate, System.IO.FileAccess.Write))
            {
                using (System.IO.StreamWriter streamWriter = new System.IO.StreamWriter(fileStream))
                {
                    WRITE(streamWriter, "   Turbine Number(m)   ", "Turbine Location-X(m)   ", "Turbine Location-Y(m)    ", "POWER(W)");
                    for (i = 1; i <= GeneralData.N_TURB; i++)
                    {
                        WRITE(streamWriter, i, GeneralData.x_turb[i - 1], GeneralData.y_turb[i - 1], GeneralData.WPOWER[i - 1]);
                    }
                }
            }
        }
        //************************************************************************	

        //************************************************************************
        //  FUNCTION : COMPUTE AREA
        //--------------------------------------------------------------
        private static void AAREA(ref double X, ref double Y, ref double Z, ref double area)
        {
            double PP;
            PP = (X + Y + Z) * 0.5;
            area = Math.Sqrt(PP * (PP - X) * (PP - Y) * (PP - Z));
            return;
        }
        //-------------------------------------------------------------
    }
}
