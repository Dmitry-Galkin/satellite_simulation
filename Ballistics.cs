using System;

namespace Satellite
{
    /// <summary>
    /// Класс Баллистика
    /// </summary>
    public class Ballistics
    {
        const double pi = 3.1415926535897932384626433832795;  // число пи
        const double a = 6378136;  // [м], большая полуось общеземного эллипсоида
        const double mu = 398600.4418e9;  // [м^3/c^2], гравитационный параметр Земли
        const double J2 = 1082.62575e-6;  // [-], коэф-т 2-ой зональной гармоники нормального потенциала
        const double eps = 1.5 * J2 * mu * a * a;  // [м^5/с], константа, определяющая сжатие Земли 
        const double EarthRate = 7.292115e-5;  // [1/c], угл. скорость Земли

        /// <summary>
        /// Наклонение орбиты, [рад]
        /// </summary>
        public double i;
        /// <summary>
        /// Cкорость изменения наклонения орбиты, [рад/с]
        /// </summary>
        public double idot;
        /// <summary>
        /// ДВУ в ИСК, [рад]
        /// </summary>
        public double lanIF;
        /// <summary>
        /// ДВУ в ГСК, [рад]
        /// </summary>
        public double lanGF;
        /// <summary>
        /// Скорость изменения ДВУ в ИСК, [рад/с]
        /// </summary>
        public double lanIFdot;  
        /// <summary>
        /// Фокальный параметр (он же радиус при круговой орбите), [м]
        /// </summary>
        public double p;
        /// <summary>
        /// Скорость изменения фокального параметра, [м/с]
        /// </summary>
        public double pdot;
        /// <summary>
        /// Элемент вектора Лапласа: k=e*sin(w), q=e*cos(w)
        /// </summary>
        public double k, q;
        /// <summary>
        /// Скорость изменения элемента вектора Лапласа
        /// </summary>
        public double kdot, qdot;  
        /// <summary>
        /// Аргумент широты, [рад]
        /// </summary>
        public double u;  
        /// <summary>
        /// Cкорость изменения аргумента широты (орбитальная угл. скорость), [рад/c]
        /// </summary>
        public double udot;  
        /// <summary>
        /// Звездное время, [рад]
        /// </summary>
        public double st;
        /// <summary>
        /// Скорость изменения звездного времени, [рад/c]
        /// </summary>
        public double stdot;  
        /// <summary>
        /// Радиус орбиты, [м]
        /// </summary>
        public double r;  
        /// <summary>
        /// Эксцентриситет орбиты
        /// </summary>
        public double e;  
        /// <summary>
        /// Оскулирующий период обращения, [с]
        /// </summary>
        public double T;
        /// <summary>
        /// Долгота, [рад]
        /// </summary>
        public double longitude;
        /// <summary>
        /// Широта, [рад]
        /// </summary>
        public double latitude;
        /// <summary>
        /// Координаты спутника в ИСК, [м]
        /// </summary>
        public double xi, yi, zi;
        /// <summary>
        /// Координаты спутника в ГСК, [м]
        /// </summary>
        public double xg, yg, zg;
        /// <summary>
        /// Матрица перехода из ИСК в ОСК
        /// </summary>
        public double[,] i2o = new double[3, 3];
        /// <summary>
        /// Матрица перехода из ГСК в ОСК
        /// </summary>
        public double[,] g2o = new double[3, 3];
        /// <summary>
        /// Кватернион ОСК->ИСК
        /// </summary>
        public double[] qoi = new double[4];
        /// <summary>
        /// Номер витка
        /// </summary>
        public int numOrbit;

        readonly double dt;  // шаг интегрирования
        readonly bool polarCompression;  // флаг учета полярного сжатия Земли
        double u0;  // угол относительно начального положения спутника       
        double S_, T_, W_;  // проекции возмущающих ускорений на оси ОСК
        int quatCoef;  // коэффициент, чтобы не было скачков в кватернионе ОСК->ИСК 
                       // при переходе аргумента широты через восходящий узел


        /// <summary>
        /// Расчет текущего радиуса орбиты
        /// </summary>
        /// <param name="p">Фокальный параметр, [м]</param>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <returns>Радиус орбиты, [м]</returns>
        public double GetOrbitRadius(double p, double q, double k, double u)
        {
            return p / (1 + q * Math.Cos(u) + k * Math.Sin(u));
        }

        /// <summary>
        /// Расчет текущей долготы 
        /// </summary>
        /// <param name="x">Координата X в ГСК</param>
        /// <param name="y">Координата Y в ГСК</param>
        /// <returns>Долгота, [рад]</returns>
        public double GetLongitude(double x, double y)
        {
            return Math.Atan2(y, x);
        }

        /// <summary>
        /// Расчет текущей широты
        /// </summary>
        /// <param name="x">Координата X в ГСК</param>
        /// <param name="y">Координата Y в ГСК</param>
        /// <param name="z">Координата Z в ГСК</param>
        /// <returns>Широта, [рад]</returns>
        public double GetLatitude(double x, double y, double z)
        {
            return Math.Asin(z / Math.Sqrt(x * x + y * y + z * z));
        }

        /// <summary>
        /// Расчет матрицы перехода из ИСК (ГСК) в ОСК (зависит от передаваемой ДВУ)
        /// </summary>
        /// <param name="lan">Долгота восходящего узла, [рад]</param>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <returns>Матрица перехода</returns>
        public double[,] GetTransMatr(double lan, double i, double u)
        {
            double[,] A = new double[3, 3];

            double slan = Math.Sin(lan);
            double clan = Math.Cos(lan);
            double si = Math.Sin(i);
            double ci = Math.Cos(i);
            double su = Math.Sin(u);
            double cu = Math.Cos(u);

            A[0, 0] = -su * clan - cu * ci * slan;
            A[0, 1] = -su * slan + cu * ci * clan;
            A[0, 2] = cu * si;
            A[1, 0] = cu * clan - su * ci * slan;
            A[1, 1] = cu * slan + su * ci * clan;
            A[1, 2] = su * si;
            A[2, 0] = -si * slan;
            A[2, 1] = si * clan;
            A[2, 2] = -ci;

            return A;
        }

        /// <summary>
        /// Расчет координат спутника в ИСК (ГСК) (зависит от передаваемой ДВУ)
        /// </summary>
        /// <param name="lan">Долгота восходящего узла, [рад]</param>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <param name="r">Радиус орбиты, [м]</param>
        /// <returns>Координаты, [м]</returns>
        public (double, double, double) GetXYZ(double lan, double i, double u, double r)
        {
            double slan = Math.Sin(lan);
            double clan = Math.Cos(lan);
            double si = Math.Sin(i);
            double ci = Math.Cos(i);
            double su = Math.Sin(u);
            double cu = Math.Cos(u);

            double x_, y_, z_;
            x_ = r * (clan * cu - slan * su * ci);
            y_ = r * (slan * cu + clan * su * ci);
            z_ = r * su * si;

            return (x_, y_, z_);
        }

        /// <summary>
        /// Расчет оскулирующего периода
        /// </summary>
        /// <param name="p">Фокальный параметр, [м]</param>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <returns>Период, [с]</returns>
        public double GetPeriod(double p, double q, double k)
        {
            // большая полуось
            double a = p / (1 - q * q - k * k);

            return 2 * pi * Math.Sqrt(a * a * a / mu);
        }

        /// <summary>
        /// Расчет эксцентриситета
        /// </summary>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <returns>Эксцентриситет</returns>
        public double GetEccentricity(double q, double k)
        {
            return Math.Sqrt(q * q + k * k);
        }

        /// <summary>
        /// Кватернион ОСК->ИСК
        /// </summary>
        /// <param name="lan">Долгота восходящего узла в ИСК, [рад]</param>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <returns>Кватернион</returns>
        public double[] GetQuat(double lan, double i, double u)
        {
            double[] q = new double[4];

            q[0] = -Math.Sin(0.5 * i) * Math.Sin(0.25 * pi + 0.5 * (lan - u));
            q[1] = Math.Cos(0.5 * i) * Math.Sin(-0.25 * pi + 0.5 * (lan + u));
            q[2] = -Math.Cos(0.5 * i) * Math.Cos(-0.25 * pi + 0.5 * (lan + u));
            q[3] = -Math.Sin(0.5 * i) * Math.Cos(0.25 * pi + 0.5 * (lan - u));

            return q;
        }

        /// <summary>
        /// Расчет производных оскулирующих элементов
        /// </summary>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="p">Фокальный параметр, [м]</param>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <param name="r">Радиус орбиты, [м]</param>
        private void OsculElemDotCalculation(double i, double p, double q, double k, double u, double r)
        {
            double si = Math.Sin(i);
            double ci = Math.Cos(i);
            double su = Math.Sin(u);
            double cu = Math.Cos(u);

            if (!polarCompression)
            {
                // орбита без возмущений
                S_ = 0;
                T_ = 0;
                W_ = 0;
            }
            else 
            {
                // орбита с возмущениями (учет полярного сжатия Земли)
                double f = -eps / (r * r * r * r);
                S_ = f * (1 - 3 * si * si * su * su);
                T_ = 2 * f * si * si * su * cu;
                W_ = 2 * f * si * ci * su;
            }

            // производные оскулирующих элементов (правые части дифуров)
            double c1 = Math.Sqrt(p / mu);
            double c2 = r / p;
            lanIFdot = c1 * c2 * su / si * W_;
            idot = c1 * c2 * cu * W_;
            pdot = 2 * c1 * r * T_; 
            qdot = c1 * (S_ * su + ((q + cu) * c2 + cu) * T_ + c2 * k * ci / si * su * W_);
            kdot = c1 * (-S_ * cu + ((k + su) * c2 + su) * T_ - c2 * q * ci / si * su * W_);
            udot = 1.0 / (c1 * c2 * r) - c1 * c2 * su * ci / si * W_;
            stdot = EarthRate;
        }

        /// <summary>
        /// Расчет всех баллистических параметров
        /// </summary>
        /// <param name="lanIF">Долгота восходящего узла в ИСК, [рад]</param>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="p">Фокальный параметр, [м]</param>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <param name="st">Звездное время, [рад]</param>
        private void BalParCalculation(double lanIF, double i, double p, double q, double k, double u, double st)
        {
            // долгота восходящего узла в ГСК
            lanGF = lanIF - st;
            if (lanGF < 0) lanGF += 2 * pi;
            // эксцентриситет
            e = GetEccentricity(q, k);
            // период
            T = GetPeriod(p, q, k);
            // радиус орбиты
            r = GetOrbitRadius(p, q, k, u);
            // координаты спутника в ИСК
            (xi, yi, zi) = GetXYZ(lanIF, i, u, r);
            // координаты спутника в ГСК
            (xg, yg, zg) = GetXYZ(lanGF, i, u, r);
            // долгота и широта
            longitude = GetLongitude(xg, yg);
            latitude = GetLatitude(xg, yg, zg);
            // матрица перехода из ИСК в ОСК
            i2o = GetTransMatr(lanIF, i, u);
            // матрица перехода из ГСК в ОСК
            g2o = GetTransMatr(lanGF, i, u);
            // кватернион
            qoi = GetQuat(lanIF, i, u);
            for (int idx = 0; idx < 4; idx++) qoi[idx] *= quatCoef;
            // производные оскулирующих элементов
            OsculElemDotCalculation(i, p, q, k, u, r);
        }

        /// <summary>
        /// Один шаг интегрирования
        /// </summary>
        public void NextStep()
        {
            // интегрирование
            lanIF += dt * lanIFdot;
            i += dt * idot;
            p += dt * pdot;
            q += dt * qdot;
            k += dt * kdot;
            u += dt * udot;
            st += dt * stdot;

            u0 += dt * udot;

            if (u >= 2 * pi)
            {
                u -= 2 * pi;
                quatCoef *= -1;
            }

            if (u0 >= 2 * pi)
            {
                // перешли на след. виток
                u0 -= 2 * pi;
                numOrbit++;
            }

            if (st >= 2 * pi) st -= 2 * pi;

            // расчет всех баллистических параметров
            BalParCalculation(lanIF, i, p, q, k, u, st);
        }

        /// <summary>
        /// Конструктор класса Баллистика
        /// </summary>
        /// <param name="lanIF">Долгота восходящего узла в ИСК, [рад]</param>
        /// <param name="i">Наклонение, [рад]</param>
        /// <param name="p">Фокальный параметр, [м]</param>
        /// <param name="q">Элемент вектора Лапласа, e*cos(w)</param>
        /// <param name="k">Элемент вектора Лапласа, e*sin(w)</param>
        /// <param name="u">Аргумент широты, [рад]</param>
        /// <param name="st">Звездное время, [рад]</param>
        /// <param name="dt">Шаг интегрирования, [с]</param>
        /// <param name="polarCompression">Флаг учета полярного сжатия Земли</param>
        public Ballistics(double lanIF = 0,
                          double i = 0,
                          double p = 7e6,
                          double q = 0,
                          double k = 0,
                          double u = 0,
                          double st = 0,
                          double dt = 0.1,
                          bool polarCompression = false)
        {
            this.lanIF = lanIF;  // ДВУ в ИСК
            this.i = i;  // наклонение
            this.p = p;  // фокальный параметр
            this.q = q;  // элемент вектора Лапласа
            this.k = k;  // элемент вектора Лапласа
            this.u = u;  // аргумент широты
            this.st = st;  // звездное время

            this.dt = dt;  // шаг интегрирования
            this.polarCompression = polarCompression;  // флаг учета полярного сжатия Земли

            numOrbit = 1;  // 1-ый виток
            u0 = 0;  // стартовая точка
            quatCoef = 1;

            // расчет всех баллистических параметров
            BalParCalculation(lanIF, i, p, q, k, u, st);
        }
    }
}
