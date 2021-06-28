using System;
using Solver;

namespace Satellite
{
    /// <summary>
    /// Класс Динамика
    /// </summary>
    public class Dynamics
    {
        /// <summary>
        /// Последовательность поворотов
        /// </summary>
        public enum RotSeq
        {
            /// <summary>
            /// psi -> gamma -> theta
            /// </summary>
            Seq213,
            /// <summary>
            /// theta -> psi -> gamma
            /// </summary>
            Seq321,
            /// <summary>
            /// psi -> theta -> gamma
            /// </summary>
            Seq231
        }
        /// <summary>
        /// Система координат
        /// </summary>
        public enum Frame
        {
            /// <summary>
            /// Орбитальная СК
            /// </summary>
            Orbital,
            /// <summary>
            /// Инерциальная СК
            /// </summary>
            ECI
        }
        /// <summary>
        /// Тензор инерции, [кг*м^2]
        /// </summary>
        public double[,] InertiaTensor
        {
            set
            {
                I = value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return I; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Ixx
        {
            set
            {
                I[0, 0] = value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return I[0, 0]; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Iyy
        {
            set
            {
                I[1, 1] = value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return I[1, 1]; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Izz
        {
            set
            {
                I[2, 2] = value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return I[2, 2]; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Ixy
        {
            set
            {
                I[0, 1] = -value;
                I[1, 0] = -value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return -I[0, 1]; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Ixz
        {
            set
            {
                I[0, 2] = -value;
                I[2, 0] = -value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return -I[0, 2]; }
        }
        /// <summary>
        /// Момент инерции, [кг*м^2]
        /// </summary>
        public double Iyz
        {
            set
            {
                I[1, 2] = -value;
                I[2, 1] = -value;
                Iinv = GetInvMatr3x3(I);
            }
            get { return -I[1, 2]; }
        }
        /// <summary>
        /// Кватернион ИСК->ССК
        /// </summary>
        public double[] qi = new double[4];
        /// <summary>
        /// Кватернион ОСК->ССК
        /// </summary>
        public double[] qo = new double[4];
        /// <summary>
        /// Угловая скорость вращения ССК относительно ИСК, [рад/с]
        /// </summary>
        public double[] wi = new double[3];
        /// <summary>
        /// Угловая скорость вращения ССК относительно ОСК, [рад/с]
        /// </summary>
        public double[] wo = new double[3];
        /// <summary>
        /// Углы ориентации ССК относительно ИСК, [рад]
        /// </summary>
        public double[] angleIF = new double[3];
        /// <summary>
        /// Углы ориентации ССК относительно ОСК, [рад]
        /// </summary>
        public double[] angleOF = new double[3];
        /// <summary>
        /// Матрица перехода из ОСК в ССК
        /// </summary>
        public double[,] o2b = new double[3, 3];
        /// <summary>
        /// Матрица перехода из ИСК в ССК
        /// </summary>
        public double[,] i2b = new double[3, 3];
        /// <summary>
        /// Суммарный кинетический момент системы, [Н*м*с]
        /// </summary>
        public double[] angMomentum = new double[3];

        double[,] I = new double[3, 3];  // тензор инерции
        double[,] Iinv = new double[3, 3];  // обращенный тензор инерции
        double[] Hrp = new double[3];  // кинетический момент вращающихся частей (маховики или гиродины, например)
        double[] Ts = new double[3];  // суммарный момент внешних сил, действующий на КА
        RotSeq rotSeq;  // последовательность поворотов
        double[] w = new double[3];  // абсолютная угловая скорость (интегрируемый параметр)
        double[] q = new double[4];  // кватернион ИСК->ССК (интегрируемый параметр)
        ODE ODESolver;  // решатель
        ODE.Method ODESolverMethod;  // метод интегрирования
        double dt;  // шаг интегрирования


        /// <summary>
        /// Обращение матрицы размера 3х3
        /// </summary>
        /// <param name="A">Обращаемая матрица</param>
        /// <returns>Обратная матрица</returns>
        private double[,] GetInvMatr3x3(double[,] A)
        {
            // обратная матрица
            double[,] B = new double[3, 3];

            // определитель матрицы
            double det = A[0, 0] * (A[1, 1] * A[2, 2] - A[2, 1] * A[1, 2]) +
                         A[0, 1] * (A[2, 0] * A[1, 2] - A[1, 0] * A[2, 2]) +
                         A[0, 2] * (A[1, 0] * A[2, 1] - A[2, 0] * A[1, 1]);
            // обращение матрицы
            B[0, 0] = (A[1, 1] * A[2, 2] - A[2, 1] * A[1, 2]) / det;
            B[1, 0] = (A[2, 0] * A[1, 2] - A[1, 0] * A[2, 2]) / det;
            B[2, 0] = (A[1, 0] * A[2, 1] - A[2, 0] * A[1, 1]) / det;

            B[0, 1] = (A[2, 1] * A[0, 2] - A[0, 1] * A[2, 2]) / det;
            B[1, 1] = (A[0, 0] * A[2, 2] - A[2, 0] * A[0, 2]) / det;
            B[2, 1] = (A[2, 0] * A[0, 1] - A[0, 0] * A[2, 1]) / det;

            B[0, 2] = (A[0, 1] * A[1, 2] - A[1, 1] * A[0, 2]) / det;
            B[1, 2] = (A[1, 0] * A[0, 2] - A[0, 0] * A[1, 2]) / det;
            B[2, 2] = (A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1]) / det;

            return B;
        }

        /// <summary>
        /// Перемножение кватернионов
        /// </summary>
        /// <param name="q1">Первый кватернион</param>
        /// <param name="q2">Второй кватернион</param>
        /// <returns>Результат перемножения</returns>
        private double[] QuatMultiplication(double[] q1, double[] q2)
        {
            double[] q_ = new double[4];

            q_[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
            q_[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
            q_[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
            q_[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

            return q_;
        }

        /// <summary>
        /// Матрица поворота от опорной СК (определяется кватернионом) к связанной СК
        /// </summary>
        /// <param name="q">Кватернион поворота от опорной СК к связанной СК</param>
        /// <returns>Матрица поворота</returns>
        public double[,] GetTransMatr(double[] q)
        {
            double[,] A = new double[3, 3];

            A[0, 0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
            A[0, 1] = 2 * (q[1] * q[2] + q[0] * q[3]);
            A[0, 2] = 2 * (q[1] * q[3] - q[0] * q[2]);
            A[1, 0] = 2 * (q[1] * q[2] - q[0] * q[3]);
            A[1, 1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
            A[1, 2] = 2 * (q[2] * q[3] + q[0] * q[1]);
            A[2, 0] = 2 * (q[1] * q[3] + q[0] * q[2]);
            A[2, 1] = 2 * (q[2] * q[3] - q[0] * q[1]);
            A[2, 2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];

            return A;
        }

        /// <summary>
        /// Вычисление углов по матрице поворота
        /// </summary>
        /// <param name="A">Матрица поворота из опорной СК в связанную СК</param>
        /// <param name="rotSeq">Последовательность поворотов</param>
        /// <returns>Углы, [рад]</returns>
        public double[] GetAngle(double[,] A, RotSeq rotSeq)
        {
            double[] angle_ = new double[3];

            // последовательность поворотов: PSI -> GAMMA -> THETA
            if (rotSeq == RotSeq.Seq213)
            {
                angle_[0] = -Math.Asin(A[2, 1]);
                angle_[1] = Math.Atan2(A[2, 0], A[2, 2]);
                angle_[2] = Math.Atan2(A[0, 1], A[1, 1]);
            }
            // последовательность поворотов: THETA -> PSI -> GAMMA
            else if (rotSeq == RotSeq.Seq321)
            {
                angle_[0] = Math.Atan2(A[1, 2], A[2, 2]);
                angle_[1] = -Math.Asin(A[0, 2]);
                angle_[2] = Math.Atan2(A[0, 1], A[0, 0]);
            }
            // последовательность поворотов: PSI -> THETA -> GAMMA
            else if (rotSeq == RotSeq.Seq231)
            {
                angle_[0] = Math.Atan2(-A[2, 1], A[1, 1]);
                angle_[1] = Math.Atan2(-A[0, 2], A[0, 0]);
                angle_[2] = Math.Asin(A[0, 1]);
            }

            return angle_;
        }

        /// <summary>
        /// Вычисление компонент кватерниона по углам
        /// </summary>
        /// <param name="angle">Углы, [рад]</param>
        /// <param name="rotSeq">Последовательность поворотов</param>
        /// <returns>Кватернион</returns>
        public double[] AngleToQuat(double[] angle, RotSeq rotSeq)
        {
            double[] q_ = new double[4];

            // последовательность поворотов: PSI -> GAMMA -> THETA
            if (rotSeq == RotSeq.Seq213)
            {
                q_[0] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) -
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[1] = Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) -
                        Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[2] = Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2) +
                        Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[3] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2) +
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2);
            }
            // последовательность поворотов: THETA -> PSI -> GAMMA
            else if (rotSeq == RotSeq.Seq321)
            {
                q_[0] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) +
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[1] = Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) -
                        Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[2] = Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2) +
                        Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[3] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2) -
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2);
            }
            // последовательность поворотов: PSI -> THETA -> GAMMA
            else if (rotSeq == RotSeq.Seq231)
            {
                q_[0] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) -
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[1] = Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Cos(angle[2] / 2) +
                        Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[2] = Math.Cos(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2) +
                        Math.Sin(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2);
                q_[3] = Math.Cos(angle[0] / 2) * Math.Cos(angle[1] / 2) * Math.Sin(angle[2] / 2) -
                        Math.Sin(angle[0] / 2) * Math.Sin(angle[1] / 2) * Math.Cos(angle[2] / 2);
            }

            return q_;
        }

        /// <summary>
        /// Расчет кинетического момента системы, [Н*м*с]
        /// </summary>
        /// <param name="w">Угловые скорости, [рад/с]</param>
        /// <param name="Hrp">Кинетический момент вращающихся частей, [Н*м*с]</param>
        /// <returns>Кинетический момент системы, [Н*м*с]</returns>
        public double[] GetAngularMomentum(double[] w, double[] Hrp = null)
        {
            double[] H_ = new double[3];

            this.Hrp = Hrp ?? new double[3] { 0, 0, 0 };
            H_[0] = I[0, 0] * w[0] + I[0, 1] * w[1] + I[0, 2] * w[2] + Hrp[0];
            H_[1] = I[1, 0] * w[0] + I[1, 1] * w[1] + I[1, 2] * w[2] + Hrp[1];
            H_[2] = I[2, 0] * w[0] + I[2, 1] * w[1] + I[2, 2] * w[2] + Hrp[2];

            return H_;
        }

        /// <summary>
        /// Правые части дифуров
        /// </summary>
        /// <param name="x">Вектор-состояния</param>
        /// <returns>Рассчитанные правые части по вектору-состояния</returns>
        private double[] RightPart(double[] x)
        {
            // абсолютные угловые скорости
            double[] w = new double[3] { x[0], x[1], x[2] };
            // кватернион ИСК->ССК
            double[] q = new double[4] { x[3], x[4], x[5], x[6] };
            
            // суммарный кинетический момент системы
            double[] Hs = GetAngularMomentum(w, Hrp);
            // промежуточная переменная для динамических ур-ний Эйлера
            double[] t = new double[3];
            t[0] = Ts[0] + Hs[1] * w[2] - Hs[2] * w[1];
            t[1] = Ts[1] + Hs[2] * w[0] - Hs[0] * w[2];
            t[2] = Ts[2] + Hs[0] * w[1] - Hs[1] * w[0];
            // вектор правых частей
            double[] rp = new double[7];
            rp[0] = Iinv[0, 0] * t[0] + Iinv[0, 1] * t[1] + Iinv[0, 2] * t[2];
            rp[1] = Iinv[1, 0] * t[0] + Iinv[1, 1] * t[1] + Iinv[1, 2] * t[2];
            rp[2] = Iinv[2, 0] * t[0] + Iinv[2, 1] * t[1] + Iinv[2, 2] * t[2];
            rp[3] = -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2]);
            rp[4] = 0.5 * (q[2] * w[2] - q[3] * w[1] + q[0] * w[0]);
            rp[5] = 0.5 * (q[3] * w[0] - q[1] * w[2] + q[0] * w[1]);
            rp[6] = 0.5 * (q[1] * w[1] - q[2] * w[0] + q[0] * w[2]);

            return rp;
        }

        /// <summary>
        /// Расчет всех параметров ориентации
        /// </summary>
        /// <param name="qoi">Кватернион ОСК->ИСК</param>
        /// <param name="orbitalAngularRate">Орбитальная угловая скорость, [рад/с]</param>
        /// <param name="rotSeq">Последовательность поворотов</param>
        private void DynParCalculation(double[] qoi, double orbitalAngularRate, RotSeq rotSeq)
        {
            // абсолютная угловая скорость
            wi = w;
            // кватернион ИСК->ССК
            qi = q;
            // кватернион ОСК->ССК
            qoi = qoi ?? new double[4] { 1, 0, 0, 0 };
            qo = QuatMultiplication(qoi, qi);
            // матрица поворота ИСК->ССК
            i2b = GetTransMatr(qi);
            // матрица поворота ОСК->ССК
            o2b = GetTransMatr(qo);
            // относительная угловая скорость КА (ОСК)
            wo[0] = wi[0] + o2b[0, 2] * orbitalAngularRate;
            wo[1] = wi[1] + o2b[1, 2] * orbitalAngularRate;
            wo[2] = wi[2] + o2b[2, 2] * orbitalAngularRate;
            // углы КА относительно ОСК
            angleIF = GetAngle(i2b, rotSeq);
            // углы КА относительно ОСК
            angleOF = GetAngle(o2b, rotSeq);
            // кинетический момент системы
            angMomentum = GetAngularMomentum(wi, Hrp);
        }

        /// <summary>
        /// Подготовка начальных условий по углам и угловым скоростям
        /// </summary>
        /// <param name="angle">Начальные углы, [рад]</param>
        /// <param name="w">Начальные угловые скорости, [рад/с]</param>
        /// <param name="qoi">Кватернион ОСК->ИСК</param>
        /// <param name="orbitalAngularRate">Орбитальная угловая скорость, [рад/с]</param>
        /// <param name="icFrame">Система координат, относиьтельно которой заданы НУ</param>
        private void InitCondPreparation(double[] angle, double[] w, double[] qoi, 
                                         double orbitalAngularRate, Frame icFrame)
        {
            // если НУ заданы относительно ОСК
            if (icFrame == Frame.Orbital)
            {
                // кватернион ОСК->ССК
                qo = AngleToQuat(angle, rotSeq);
                // матрица поворота ОСК -> ССК
                o2b = GetTransMatr(qo);
                // абсолютная угловая скорость КА
                this.w[0] = w[0] - o2b[0, 2] * orbitalAngularRate;
                this.w[1] = w[1] - o2b[1, 2] * orbitalAngularRate;
                this.w[2] = w[2] - o2b[2, 2] * orbitalAngularRate;
                // нахождение кватерниона ИСК->ССК
                q[0] = qoi[0];
                q[1] = -qoi[1];
                q[2] = -qoi[2];
                q[3] = -qoi[3];
                this.q = QuatMultiplication(q, qo);
            }
            // если НУ заданы относительно ИСК
            else if (icFrame == Frame.ECI)
            {
                // начальные угловая скорость
                this.w = w;
                // кватернион ИСК->ССК
                this.q = AngleToQuat(angle, rotSeq);
            }
        }

        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        /// <param name="Ts">Суммарный момент внешних сил, [Н*м]/param>
        /// <param name="Hrp">Кинетический момент вращающихся частей, [Н*м*с]</param>
        /// <param name="qoi">Кватернион ОСК->ИСК</param>
        /// <param name="orbitalAngularRate">Орбитальная угловая скорость</param>
        /// <param name="dt">Шаг интегрирования, [с]</param>
        public void NextStep(double[] Ts,
                             double[] Hrp = null,
                             double[] qoi = null,
                             double orbitalAngularRate = 0,
                             double dt = -1)
        {
            // суммарный момент внешних сил, действующий на КА
            this.Ts = Ts;
            // кинетический момент вращающихся частей (маховики или гиродины, например)
            this.Hrp = Hrp ?? new double[3] { 0, 0, 0 };
            // НУ
            double[] x0 = new double[7] { w[0], w[1], w[2], q[0], q[1], q[2], q[3] };
            // решение системы дифуров
            if (dt < 0) dt = this.dt;
            double[] y = ODESolver.NextStep(x0, dt);
            w[0] = y[0];
            w[1] = y[1];
            w[2] = y[2];
            q[0] = y[3];
            q[1] = y[4];
            q[2] = y[5];
            q[3] = y[6];
            // нормировка кватерниона
            double L = Math.Sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            q[0] /= L;
            q[1] /= L;
            q[2] /= L;
            q[3] /= L;
            // расчет всех параметров ориентации
            DynParCalculation(qoi, orbitalAngularRate, rotSeq);
        }

        /// <summary>
        /// Конструктор класс Динамика
        /// </summary>
        /// <param name="Jxx">Момент инерции, [кг*м^2]</param>
        /// <param name="Jyy">Момент инерции, [кг*м^2]</param>
        /// <param name="Jzz">Момент инерции, [кг*м^2]</param>
        /// <param name="Jxy">Момент инерции, [кг*м^2]</param>
        /// <param name="Jxz">Момент инерции, [кг*м^2]</param>
        /// <param name="Jyz">Момент инерции, [кг*м^2]</param>
        /// <param name="angle">Начальные углы, [рад]</param>
        /// <param name="w">Начальные угловые скорости, [рад/с]</param>
        /// <param name="icFrame">СК, относительно которой заданы НУ</param>
        /// <param name="rotSeq">Последовательность поворотов</param>
        /// <param name="qoi">Кватернион ОСК->ИСК</param>
        /// <param name="orbitalAngularRate">Орбитальная угловая скорость, [рад/с]</param>
        /// <param name="ODESolverMethod">Метод решения дифуров</param>
        /// <param name="dt">Шаг интегрирования, [с]</param>
        public Dynamics(double Jxx = 1, double Jyy = 1, double Jzz = 1,
                        double Jxy = 0, double Jxz = 0, double Jyz = 0,
                        double[] angle = null, double[] w = null, double[] Hrp = null,
                        Frame icFrame = Frame.Orbital, RotSeq rotSeq = RotSeq.Seq231,
                        double[] qoi = null,
                        double orbitalAngularRate = 0,
                        ODE.Method ODESolverMethod = ODE.Method.RK4,
                        double dt = 0.1)
        {
            // тензор инерции
            Ixx = Jxx;
            Iyy = Jyy;
            Izz = Jzz;
            Ixy = Jxy;
            Ixz = Jxz;
            Iyz = Jyz;
            // последовательность поворотов
            this.rotSeq = rotSeq;
            // метод интегрирования
            this.ODESolverMethod = ODESolverMethod;
            // шаг интегрирования
            this.dt = dt;

            // подготовка начальных условий по углам и угловым скоростям
            angle = angle ?? new double[3] { 0, 0, 0 };
            w = w ?? new double[3] { 0, 0, 0 };
            InitCondPreparation(angle, w, qoi, orbitalAngularRate, icFrame);

            // решатель дифуров динамики
            ODESolver = new ODE(ODESolverMethod, RightPart, dt);

            // расчет всех параметров ориентации
            this.Hrp = Hrp ?? new double[3] { 0, 0, 0 };
            DynParCalculation(qoi, orbitalAngularRate, rotSeq);
        }
    }
}
