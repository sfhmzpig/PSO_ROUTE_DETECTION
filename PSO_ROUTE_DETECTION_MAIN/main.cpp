#include <iostream>
#include <ctime>
#include <fstream>
#include <iterator>
#include <cmath>
const int X_max = 100;
const int X_min = 0;
const int Y_max = 100;
const int Y_min = 0;
const int Route_number = 1;
const int Route_point = 101;
const int Dimension = 2;
const double PI = 3.1415926;

class R
{
private:
    double res;
    double NUMDEV;
    double NUMMIN;
    double NUMMAX;
    double calRes()
    {
        res = rand()/(double)RAND_MAX;
    }

public:
    R()
    {
        srand((unsigned)time(0));
    }
    double get_res(int min,int max)
    {
        calRes();
        return res * (max-min) + min;
    }

};

class Route
{
private:
    int route_number;
    int route_point_number;
    double ***route;

public:
    Route(std::string map_file_name,int route_num,int route_point_num)
    {
        route_number = route_num;
        route_point_number = route_point_num;

        route = new double**[route_number];
        for (int i = 0; i != route_number; ++i)
        {
            route[i] = new double*[route_point_number];
            for (int j = 0 ; j != route_point_number; ++j)
            {
                route[i][j] = new double[2];
            }
        }

        std::ifstream map_stream(map_file_name);
        std::istream_iterator<double> route_of_map_ite(map_stream);
        for (int i = 0; i != route_number; ++i)
        {
            for (int j = 0; j != route_point_number; ++j)
            {
                for (int k = 0; k != Dimension; ++k)
                {
                    route[i][j][k] = *route_of_map_ite++;
                }
            }
        }
        std::cout << "ROUTE ok" << std::endl;
    }

    double *** get_route()
    {
        double ***tmp_route = route;
        return tmp_route;
    }
};


class Fit
{
private:
    int route_number;
    int point_number;
    bool **cover_bit;
    double cover_ratio;
public:
    Fit(int route_num,int point_num)
    {
        route_number = route_num;
        point_number = point_num;
        cover_bit = new bool*[route_number];
        for (int i = 0; i != route_number; ++i)
        {
            cover_bit[i] = new bool[point_number];
        }
        for (int i = 0; i != route_number; ++i)
        {
            for (int j = 0; j != point_number; ++j)
                cover_bit[i][j] = false;
        }
        std::cout << "FIT ok" << std::endl;

    }
    void set_cover_bit(int route_index,int point_index)
    {
        cover_bit[route_index][point_index] = true;
    }
    void cal_ratio()
    {
        int cover_num = 0;
        for (int i = 0; i != route_number; ++i)
        {
            for (int j = 0; j != point_number; ++j)
            {
                if (cover_bit[i][j])
                    ++cover_num;
            }
        }
        cover_ratio = (double)cover_num/(route_number*point_number);
    }
    double get_cover_ratio()
    {
        return cover_ratio;
    }
    void reset_cover_bit()
    {
        for (int i = 0; i != route_number; ++i)
        {
            for (int j = 0; j != point_number; ++j)
                cover_bit[i][j] = false;
        }
    }
};

class PSO
{
private:
    double ***person;    //个体位置
    double *every_person_fit;
    double ***person_v;  //个体速度
    double *p_best; //个体最优
    double ***p_best_location;
    double *g_best;  //群体最优
    double **g_best_location;
    double v_max;
    double v_min;
    double w_min;
    double w_max;
    double w;
    double c1;
    double c2;
    int max_person_number;
    int max_intration;
    int index_interation;
    int D;
    int radius;
    R r = R();
    Route route = Route("map.txt",Route_number,Route_point);
    Fit fit = Fit(Route_number,Route_point);
public:
    PSO()
    {
        radius = 2;
        v_max = 100;
        v_min = -100;
        w_min = 0.4;
        w_max = 0.9;
        max_person_number = 100;
        max_intration = 1000;
        D = X_max+1;
        c1 = 2.0;
        c2 = 2.0;
        std::cout << "PSO ok" << std::endl;

    }

    void initialize()
    {
        person = new double**[max_person_number];
        for (int i = 0; i != max_person_number; ++i)
        {
            person[i] = new double*[D];
            for (int j = 0; j != D; ++j)    //二维空间
                person[i][j] = new double[Dimension];
        }
        person_v = new double**[max_person_number];
        for (int i = 0; i != max_person_number; ++i)
        {
            person_v[i] = new double*[D];
            for (int j = 0; j != D; ++j)    //二维空间
                person_v[i][j] = new double[Dimension];
        }

        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                for (int k = 0 ; k != Dimension; ++k)
                {
                    person[i][j][k] = r.get_res(X_min,X_max);
                    //std::cout << person[i][j][k] << std::endl;
                }
            }
        }

        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                for (int k = 0 ; k != Dimension; ++k)
                {
                    person_v[i][j][k] = r.get_res((int)v_min,(int)v_max);
                    //std::cout << person_v[i][j][k] << std::endl;
                }
            }
        }

        every_person_fit = new double[max_person_number];
        every_person_fit = calEveryPersonFit();
        p_best = new double[max_person_number];
        for (int i = 0; i != max_person_number; ++i)
        {
            p_best[i] = every_person_fit[i];
            //std::cout << p_best[i] << std::endl;
        }
        g_best = new double[2];
        g_best = select_best();
        g_best_location = new double*[D];
        for (int i = 0; i != D; ++i)
        {
            g_best_location[i] = new double[Dimension];
        }
        for (int i = 0; i != D; ++i)
        {
            for (int j = 0; j != Dimension; ++j)
            {
                g_best_location[i][j] = person[(int)g_best[1]][i][j];
            }
            //std::cout << g_best_location[i][0] << '\t' << g_best_location[i][1] << std::endl;
        }

        p_best_location = new double **[max_person_number];
        for (int i = 0; i != max_person_number; ++i)
        {
            p_best_location[i] = new double *[D];
            for (int j = 0; j != D; ++j)
            {
                p_best_location[i][j] = new double[Dimension];
            }
        }
        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                for (int k = 0; k != Dimension; ++k)
                {
                    p_best_location[i][j][k] = person[i][j][k];
                }
                //std::cout << p_best_location[i][j][0] << '\t' << p_best_location[i][j][1] << std::endl;
            }
        }
        std::cout << "PSO::initialize ok" << std::endl;

    }

    void one_evolution(int count)
    {
        w = set_w(count);
        //std::cout << w << std::endl;
        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                for (int k = 0; k != Dimension; ++k)
                {
                    //std::cout << person[i][j][k] << std::endl;
                    person_v[i][j][k] = w * person_v[i][j][k] + c1 * r.get_res(0,1) * (g_best_location[j][k] - person[i][j][k]) + c2 * r.get_res(0,1) * (p_best_location[i][j][k] - person[i][j][k]);
                    if (person_v[i][j][k] > v_max)
                        person_v[i][j][k] = v_max;
                    if (person_v[i][j][k] < v_min)
                        person_v[i][j][k] = v_min;
                    //std::cout << person_v[i][j][k] << std::endl;
                }
            }
        }
        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                for (int k = 0; k != Dimension; ++k)
                {
                    person[i][j][k] += person_v[i][j][k];
                    if (person[i][j][k] > X_max)
                        person[i][j][k] = X_max;
                    if (person[i][j][k] < X_min)
                        person[i][j][k] = X_min;
                    //std::cout << person[i][j][k] << std::endl;
                }
            }
        }
        std::cout << "evolution::location ok" << std::endl;

        every_person_fit = calEveryPersonFit();
        /*
        for (int i = 0 ; i != max_person_number; ++i)
        {
            std::cout << every_person_fit[i] << std::endl;
        }
         */
        std::cout << "evolution::location ok" << std::endl;
        for (int i = 0; i != max_person_number; ++i)
        {
            if (every_person_fit[i] > p_best[i])
            {
                p_best[i] = every_person_fit[i];
                for (int j = 0; j != D; ++j)
                {
                    for (int k = 0 ; k != Dimension; ++k)
                    {
                        p_best_location[i][j][k] = person[i][j][k];
                    }
                }
            }
        }
        std::cout << "evolution::p ok" << std::endl;

        g_best = select_best();
        std::cout << g_best[0] << '\t' << g_best[1] << std::endl;

        for (int i = 0; i != D; ++i)
        {
            for (int j = 0; j != Dimension; ++j)
            {
                g_best_location[i][j] = person[(int)g_best[1]][i][j];
            }
            //std::cout << g_best_location[i][0] << '\t' << g_best_location[i][1] << std::endl;
        }

        std::cout << "evolution::g ok" << std::endl;

        std::cout << "PSO::evolution " << count << " ok" << std::endl;


    }

    double * select_best()
    {
        double best_fit = p_best[0];
        int select = 0;
        for (int i = 1; i != max_person_number; ++i)
        {
            if (p_best[i] >= best_fit)
            {
                best_fit = p_best[i];
                select = i;
            }
        }
        g_best[0] = best_fit;
        g_best[1] = select;
        return g_best;
    }

    double set_w(int count)
    {
        double w;
        w = (w_max - w_min)*(max_intration - count)/max_intration + w_min;
        return w;
    }

    double * calEveryPersonFit()
    {
        double *tmp_every_person_fit = new double[max_person_number];
        /*
        double **tmp_every_person_loc = new double*[max_person_number];
        for (int i = 0 ; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                tmp_every_person_loc[i] = new double[2];
            }
        }
        */
        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                double tmp_x = person[i][j][0];
                double tmp_y = person[i][j][1];
                double *** tmp_route = route.get_route();
                for (int k = 0; k != Route_number; ++k)
                {
                    for (int ii = 0; ii != Route_point; ++ii)
                    {
                        double tmp_distance_of_person_and_route_point = sqrt((tmp_x - tmp_route[k][ii][0])*(tmp_x - tmp_route[k][ii][0]) + (tmp_y - tmp_route[k][ii][1])*(tmp_y - tmp_route[k][ii][1]));
                        if (tmp_distance_of_person_and_route_point <= radius)
                        {
                            fit.set_cover_bit(k,ii);
                        }
                    }
                }
            }
            fit.cal_ratio();
            tmp_every_person_fit[i] = fit.get_cover_ratio();
            fit.reset_cover_bit();
        }
        /*
        double *** tmp_eve_per_fit_and_loc;
        double ** tmp_eve_per_fit_package = new double*;
        tmp_eve_per_fit_package[0] = tmp_every_person_fit;
        tmp_eve_per_fit_and_loc[0] = tmp_eve_per_fit_package;
        tmp_eve_per_fit_and_loc[1] = person
         */
        return tmp_every_person_fit;
    }

    void cout_person(std::string result_file_name)
    {
        std::ofstream result_stream(result_file_name);

        for (int i = 0; i != max_person_number; ++i)
        {
            for (int j = 0; j != D; ++j)
            {
                result_stream << person[i][j][0] << '\t' << person[i][j][1] << std::endl;
            }
        }
    }

    void PSO_begin()
    {
        initialize();
        for (int i = 0; i != max_intration; ++i)
        {
            one_evolution(i);
        }
        cout_person("result.txt");
    }
};




int main() {
    PSO start_pso = PSO();
    start_pso.PSO_begin();
    return 0;
}