namespace Solver
{
    public class ODE
    {
        /// <summary>
        /// Метод численного интегрирования
        /// </summary>
        public enum Method
        {
            /// <summary>
            /// Метод Эйлера 1-го поядка
            /// </summary>
            Euler,
            /// <summary>
            /// Метод Рунге-Кутты 2-го порядка (модиф. метод Эйлера)
            /// </summary>
            RK2,
            /// <summary>
            /// Метод Рунге-Кутты 4-го порядка
            /// </summary>
            RK4
        }

        public delegate double[] RightPart(double[] x);

        private readonly Method method;  // метод интегрирования
        private readonly RightPart rp;  // правые части
        private readonly double dt;  // шаг интегрирования


        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        /// <param name="x0">Начальные условия</param>
        /// <returns>Решение в конце интервала</returns>
        public double[] NextStep(double[] x0)
        {
            return NextStep_(method, rp, x0, dt);
        }

        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        public double[] NextStep(double[] x0, double dt)
        {
            return NextStep_(method, rp, x0, dt);
        }

        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        /// <param name="method">Метод интегрирования</param>
        /// <param name="rp">Правые части</param>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        public double[] NextStep(Method method, RightPart rp, double[] x0, double dt)
        {
            return NextStep_(method, rp, x0, dt);
        }

        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        /// <param name="method">Метод интегрирования</param>
        /// <param name="rp">Правые части</param>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        private double[] NextStep_(Method method, RightPart rp, double[] x0, double dt)
        {
            double[] v = new double[x0.Length];

            switch (method)
            {
                case Method.Euler:
                    v = EulerSolver(rp, x0, dt);
                    break;
                case Method.RK2:
                    v = RK2Solver(rp, x0, dt);
                    break;
                case Method.RK4:
                    v = RK4Solver(rp, x0, dt);
                    break;
            }

            return v;
        }

        /// <summary>
        /// Метод Эйлера
        /// </summary>
        /// <param name="rp">Правые части</param>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        private double[] EulerSolver(RightPart rp, double[] x0, double dt)
        {
            double[] v = new double[x0.Length];  // выход

            double[] k = rp(x0);  // правые части
            for (int i = 0; i < x0.Length; i++)
            {
                v[i] = x0[i] + dt * k[i];
            }

            return v;
        }

        /// <summary>
        /// Метод Рунге-Кутты 2-го порядка (модиф. метод Эйлера)
        /// </summary>
        /// <param name="rp">Правые части</param>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        private double[] RK2Solver(RightPart rp, double[] x0, double dt)
        {
            double[] v = new double[x0.Length];  // выход

            double[] midpoint = new double[x0.Length];
            // ПРОГНОЗ
            double[] k1 = rp(x0);
            for (int i = 0; i < x0.Length; i++)
            {
                midpoint[i] = x0[i] + dt * k1[i];
            }
            // КОРРЕКЦИЯ
            double[] k2 = rp(midpoint);
            for (int i = 0; i < x0.Length; i++)
            {
                v[i] = x0[i] + dt * (k1[i] + k2[i]) / 2;
            }

            return v;
        }

        /// <summary>
        /// Метод Рунге-Кутты 4-го порядка
        /// </summary>
        /// <param name="rp">Правые части</param>
        /// <param name="x0">Начальные условия</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        /// <returns>Решение в конце интервала</returns>
        private double[] RK4Solver(RightPart rp, double[] x0, double dt)
        {
            double[] v = new double[x0.Length];  // выход

            double[] midpoint = new double[x0.Length];
            // k1 - наклон в начале интервала
            double[] k1 = rp(x0);
            for (int i = 0; i < x0.Length; i++)
            {
                k1[i] *= dt;
                midpoint[i] = x0[i] + 0.5 * k1[i];
            }
            // k2 - наклон в промежуточной точке интервала
            double[] k2 = rp(midpoint);
            for (int i = 0; i < x0.Length; i++)
            {
                k2[i] *= dt;
                midpoint[i] = x0[i] + 0.5 * k2[i];
            }
            // k3 - наклон в промежуточной точке интервала
            double[] k3 = rp(midpoint);
            for (int i = 0; i < x0.Length; i++)
            {
                k3[i] *= dt;
                midpoint[i] = x0[i] + k3[i];
            }
            // k4 - наклон в конце интервала
            double[] k4 = rp(midpoint);

            for (int i = 0; i < x0.Length; i++)
            {
                k4[i] *= dt;
                v[i] = x0[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
            }

            return v;
        }

        /// <summary>
        /// ОДУ
        /// </summary>
        /// <param name="method">Метод интегрирования</param>
        /// <param name="rp">Правые части</param>
        /// <param name="dt">Шаг интегрирования, [c]</param>
        public ODE(Method method = Method.Euler, RightPart rp = null, double dt = 0.1)
        {
            // метод интегрирования
            this.method = method;
            // правые части
            this.rp = rp;
            // шаг интегрирования
            this.dt = dt;
        }
    }
}
