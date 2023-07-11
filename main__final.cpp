/*********************Alekseev A. V.*******************************
************************************************************/
#include <iostream>
#include <iomanip>
#include <math.h>


using namespace std;

//структура для возврата координат точки из объекта Vector3D
struct Coords { double x, y, z; };
//структура для возврата точки пересечения отрезков из функции intersect
struct Intersect{ double x,y,z; int flag; };

class Vector3D //класс вектора, обозначающего конец отрезка
{
    public:
        Vector3D(double x_in=0, double y_in=0, double z_in=0)
        {
            x=x_in;
            y=y_in;
            z=z_in;
        }
        Coords get_point_coords(); //получить координаты точки конца вектора, хранящиеся в объекте
        const Vector3D &operator=(const Vector3D &v3d); //перегрузка оператора присваивания

    protected: //чтобы можно было реализовать наследование
        double x;
        double y;
        double z;

};
Coords Vector3D::get_point_coords()
{
    Coords  current_coords;
            current_coords.x=x;
            current_coords.y=y;
            current_coords.z=z;
    return current_coords;
}
const Vector3D &Vector3D::operator=(const Vector3D &v3d)
{
   if(this!=&v3d){
        x = v3d.x;
        y = v3d.y;
        z = v3d.z;
   }
   return *this;
}



class Segment3D  //класс отрезка построенного по двум векторам
{
    public:
        Segment3D(Vector3D v_start_in, Vector3D v_end_in){ v_start=v_start_in; v_end=v_end_in; }
        Coords get_point_coords_defined_by_parameter(double parameter); //получить координаты точки, пренадлежащей отрезку, и заданной параметром в параметрических уравнениях прямой
        Coords get_start_coords(); //получить координаты начала отрезка
        Coords get_end_coords();   //получить координаты конца отрезка
    protected:
        Vector3D v_start;
        Vector3D v_end;
};
Coords Segment3D::get_point_coords_defined_by_parameter(double parameter)
{
    Coords point_coords, seg_start, seg_end;

    seg_start=get_start_coords();
    seg_end=get_end_coords();

    point_coords.x=seg_start.x + parameter*(seg_end.x - seg_start.x);
    point_coords.y=seg_start.y + parameter*(seg_end.y - seg_start.y);
    point_coords.z=seg_start.z + parameter*(seg_end.z - seg_start.z);

    return point_coords;
}

Coords Segment3D::get_start_coords()
{
   return v_start.get_point_coords();
}

Coords Segment3D::get_end_coords()
{
   return v_end.get_point_coords();
}



//  *****************************************************************************************************************
//функция получает два объекта отрезков и возвращает структуру Intersect, содержащую
//координаты точки пересечения отрезков и флаг взаимного положения отрезков
// если отрезки пересекаются то flag=1, если совпадают flag=2, если не пересекаются flag=3, есди ошибка при расчетах flag=4
Intersect intersect(Segment3D seg1, Segment3D seg2)
{
    double system_extended_matrix [3][3];   //расширенная матрица системы уравнений для пересечения отрезков
    int r = -1;                             //ранг матрицы системы
    int r_extend = -1;                      //ранг расширенной матрицы системы
    double t,f;                             //параметры для вычисления точки пересечения отрезков по параметрическим уравнениям
    double delta_kramer, delta_t_kramer, delta_f_kramer; //определители матрицы системы уравнений для решения методом Крамера
    Coords intersect_point_seg1, intersect_point_seg2; //переменные для хранения результатов вычисления коорднат по параметрам для каждого отрезка
    Intersect result;                       //структура для возврата результата работы функции
    result.x=0;
    result.y=0;
    result.z=0;
    result.flag=0;
    double extended_matrix_determinante;    //определитель для расширенной матрицы системы
    double minors[3][3];                    //миноры для расширенной матрицы системы и для простой матрицы системы
                        //первый столбец -- миноры, одинаковые для простой матрицы и расширенной
                        //1 и 2 столбец миноры расширенной матрицы включающие свободные члены,
                        //полученные из первого минора текущей строки, заменой его столбцов на свободные члены
                        //первый столбец -- замена первого столбца минора, стоящего в 0 столбце этой строки
                        //второй столбец -- замена второго столбца минора, стоящего в 0 столбце этой строки

        //коэффициенты уравнений системы перед параметром t (относятся к первому отрезку)
    system_extended_matrix[0][0] = seg1.get_end_coords().x - seg1.get_start_coords().x;     //a00
    system_extended_matrix[1][0] = seg1.get_end_coords().y - seg1.get_start_coords().y;     //a10
    system_extended_matrix[2][0] = seg1.get_end_coords().z - seg1.get_start_coords().z;     //a20
        //коэффициенты уравнений системы перед параметром f (относятся ко второму отрезку)
    system_extended_matrix[0][1] = seg2.get_start_coords().x - seg2.get_end_coords().x;     //a01
    system_extended_matrix[1][1] = seg2.get_start_coords().y - seg2.get_end_coords().y;     //a11
    system_extended_matrix[2][1] = seg2.get_start_coords().z - seg2.get_end_coords().z;     //a21
        //коэффициенты уравнений ситемы отноятся к свободным членам
    system_extended_matrix[0][2] = seg2.get_start_coords().x - seg1.get_start_coords().x; //c1=a02
    system_extended_matrix[1][2] = seg2.get_start_coords().y - seg1.get_start_coords().y; //c2=a12
    system_extended_matrix[2][2] = seg2.get_start_coords().z - seg1.get_start_coords().z; //c3=a22

    //перебор элементов расширенной матрицы:
        //если все =0 то отрезки совпадают ---ВЫХОД
    for(int i=0; i<3; i++){
        for(int j=0; j<3;j++){
            if(system_extended_matrix[i][j] != 0 ) {  //здесь правильнее задать не ноль а минимальную погрешность, но в условиях задачи погрешность не указана
                    r_extend = 1;
                    break;
            }
        }
        if(r_extend==1) break;
    }
    if(r_extend == -1) { result.flag=2; return result; }

    //считаем определитель расширенной матрицы:
        //если не равен нолю значит ранг 3 --- система несовместна, отрезки не пересекаются  -- ВЫХОД
    extended_matrix_determinante = system_extended_matrix[0][0] * system_extended_matrix[1][1] * system_extended_matrix[2][2]
                                   + system_extended_matrix[1][0] *  system_extended_matrix[2][1] * system_extended_matrix[0][2]
                                   + system_extended_matrix[0][1] *  system_extended_matrix[1][2] * system_extended_matrix[2][0]
                                   - system_extended_matrix[0][2] *  system_extended_matrix[1][1] * system_extended_matrix[2][0]
                                   - system_extended_matrix[0][1] *  system_extended_matrix[1][0] * system_extended_matrix[2][2]
                                   - system_extended_matrix[1][2] *  system_extended_matrix[2][1] * system_extended_matrix[0][0];
    if(extended_matrix_determinante != 0) { result.flag=3; return result;  }

  //вычисляем миноры матрицы
    minors[0][0] = system_extended_matrix[0][0] * system_extended_matrix[1][1] -
                    system_extended_matrix[0][1] * system_extended_matrix[1][0];
    minors[1][0] = system_extended_matrix[1][0] * system_extended_matrix[2][1] -
                    system_extended_matrix[1][1] * system_extended_matrix[2][0];
    minors[2][0] = system_extended_matrix[0][0] * system_extended_matrix[2][1] -
                    system_extended_matrix[0][1] * system_extended_matrix[2][0];

    minors[0][1] = system_extended_matrix[0][2] * system_extended_matrix[1][1] -
                    system_extended_matrix[1][2] * system_extended_matrix[0][1];
    minors[0][2] = system_extended_matrix[0][0] * system_extended_matrix[1][2] -
                    system_extended_matrix[1][0] * system_extended_matrix[0][2];
    minors[1][1] = system_extended_matrix[1][2] * system_extended_matrix[2][1] -
                    system_extended_matrix[1][1] * system_extended_matrix[2][2];
    minors[1][2] = system_extended_matrix[1][0] * system_extended_matrix[2][2] -
                    system_extended_matrix[1][2] * system_extended_matrix[2][0];
    minors[2][1] = system_extended_matrix[0][2] * system_extended_matrix[2][1] -
                    system_extended_matrix[0][1] * system_extended_matrix[2][2];
    minors[2][2] = system_extended_matrix[0][0] * system_extended_matrix[2][2] -
                    system_extended_matrix[0][2] * system_extended_matrix[2][0];

    //перебираем миноры простой матрицы и определяем ее ранг, а заодно запоминаем миноры для решения системы методом Крамера
    for(int i=0; i<3; i++){
       if( minors[i][0] != 0) { //здесь правильнее задать не ноль а минимальную погрешность, но в условиях задачи погрешность не указана
            r=2;
            delta_kramer = minors[i][0];
            delta_t_kramer = minors[i][1];
            delta_f_kramer = minors[i][2];
            break;
       }
       else r=1;
    }

    //перебираем миноры расширенной матрицы и определяем ее ранг
    for(int i=0; i<3; i++){
        for(int j=1;j<3;j++){
            if( minors[i][j] != 0) { r_extend=2; break; } //здесь правильнее задать не ноль а минимальную погрешность, но в условиях задачи погрешность не указана
            else r_extend=1;
        }
        if(r_extend == 2) break;
    }
    r_extend = max (r, r_extend);

    //сравниваем r и r_extend если не равны то отрезки не пересекаются -- ВЫХОД
    if(r != r_extend) { result.flag=3; return result; }

    //если равны но меньше 2 --- отрезки совпадают --- ВЫХОД
    if(r_extend <2 ) { result.flag=1; return result; }

    //если равны и равны 2 то решаем систему методом Крамера используя рассчитанные раньше определители
    //если любое другое число -- ошибка функции
    if(r_extend != 2) { result.flag=4; return result; }

    t= delta_t_kramer / delta_kramer;
    f= delta_f_kramer / delta_kramer;

    //получаем для seg1 координаты точки на нем, определеннные параметром t
    intersect_point_seg1=seg1.get_point_coords_defined_by_parameter(t);

    //получаем для seg2 координаты точки на нем, определеннные параметром f
    intersect_point_seg2=seg2.get_point_coords_defined_by_parameter(f);

    //сравниваем полученные координаты точки пересечения отрезков seg1 seg 2, они должны совпадать
    if(intersect_point_seg1.x == intersect_point_seg2.x &&
       intersect_point_seg1.y == intersect_point_seg2.y &&
       intersect_point_seg1.z == intersect_point_seg2.z){
            result.flag=1;
            result.x=intersect_point_seg1.x;
            result.y=intersect_point_seg1.y;
            result.z=intersect_point_seg1.z;
            return result;
    }
    else { result.flag=4; return result; }

}


// __******************************************** драйвер для проверки функции ********************************************
int main()
{
    Vector3D
        seg1_start(1.1, 1, 0),
        seg1_end(7.1, 7, 1),
        seg2_start(1.1, 7, 0),
        seg2_end(7.1, 1, 1);


    Segment3D
        seg1(seg1_start, seg1_end),
        seg2(seg2_start, seg2_end);

    Intersect result;

        result = intersect(seg1, seg2);

        switch(result.flag){
            case 1: cout<<"Intersect coords:"<<endl<<" x="<<result.x<<" y="<<result.y<<" z="<<result.z; break;
            case 2: cout<<"segments overlaps";break;
            case 3: cout<<"segments have no points of intersects"; break;
            case 4: cout<<"function calculation error";break;

            default: cout<<"error";break;
        }


    return 0;
}


