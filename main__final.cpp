/*********************Alekseev A. V.*******************************
************************************************************/
#include <iostream>
#include <iomanip>
#include <math.h>


using namespace std;

//��������� ��� �������� ��������� ����� �� ������� Vector3D
struct Coords { double x, y, z; };
//��������� ��� �������� ����� ����������� �������� �� ������� intersect
struct Intersect{ double x,y,z; int flag; };

class Vector3D //����� �������, ������������� ����� �������
{
    public:
        Vector3D(double x_in=0, double y_in=0, double z_in=0)
        {
            x=x_in;
            y=y_in;
            z=z_in;
        }
        Coords get_point_coords(); //�������� ���������� ����� ����� �������, ���������� � �������
        const Vector3D &operator=(const Vector3D &v3d); //���������� ��������� ������������

    protected: //����� ����� ���� ����������� ������������
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



class Segment3D  //����� ������� ������������ �� ���� ��������
{
    public:
        Segment3D(Vector3D v_start_in, Vector3D v_end_in){ v_start=v_start_in; v_end=v_end_in; }
        Coords get_point_coords_defined_by_parameter(double parameter); //�������� ���������� �����, ������������� �������, � �������� ���������� � ��������������� ���������� ������
        Coords get_start_coords(); //�������� ���������� ������ �������
        Coords get_end_coords();   //�������� ���������� ����� �������
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
//������� �������� ��� ������� �������� � ���������� ��������� Intersect, ����������
//���������� ����� ����������� �������� � ���� ��������� ��������� ��������
// ���� ������� ������������ �� flag=1, ���� ��������� flag=2, ���� �� ������������ flag=3, ���� ������ ��� �������� flag=4
Intersect intersect(Segment3D seg1, Segment3D seg2)
{
    double system_extended_matrix [3][3];   //����������� ������� ������� ��������� ��� ����������� ��������
    int r = -1;                             //���� ������� �������
    int r_extend = -1;                      //���� ����������� ������� �������
    double t,f;                             //��������� ��� ���������� ����� ����������� �������� �� ��������������� ����������
    double delta_kramer, delta_t_kramer, delta_f_kramer; //������������ ������� ������� ��������� ��� ������� ������� �������
    Coords intersect_point_seg1, intersect_point_seg2; //���������� ��� �������� ����������� ���������� �������� �� ���������� ��� ������� �������
    Intersect result;                       //��������� ��� �������� ���������� ������ �������
    result.x=0;
    result.y=0;
    result.z=0;
    result.flag=0;
    double extended_matrix_determinante;    //������������ ��� ����������� ������� �������
    double minors[3][3];                    //������ ��� ����������� ������� ������� � ��� ������� ������� �������
                        //������ ������� -- ������, ���������� ��� ������� ������� � �����������
                        //1 � 2 ������� ������ ����������� ������� ���������� ��������� �����,
                        //���������� �� ������� ������ ������� ������, ������� ��� �������� �� ��������� �����
                        //������ ������� -- ������ ������� ������� ������, �������� � 0 ������� ���� ������
                        //������ ������� -- ������ ������� ������� ������, �������� � 0 ������� ���� ������

        //������������ ��������� ������� ����� ���������� t (��������� � ������� �������)
    system_extended_matrix[0][0] = seg1.get_end_coords().x - seg1.get_start_coords().x;     //a00
    system_extended_matrix[1][0] = seg1.get_end_coords().y - seg1.get_start_coords().y;     //a10
    system_extended_matrix[2][0] = seg1.get_end_coords().z - seg1.get_start_coords().z;     //a20
        //������������ ��������� ������� ����� ���������� f (��������� �� ������� �������)
    system_extended_matrix[0][1] = seg2.get_start_coords().x - seg2.get_end_coords().x;     //a01
    system_extended_matrix[1][1] = seg2.get_start_coords().y - seg2.get_end_coords().y;     //a11
    system_extended_matrix[2][1] = seg2.get_start_coords().z - seg2.get_end_coords().z;     //a21
        //������������ ��������� ������ �������� � ��������� ������
    system_extended_matrix[0][2] = seg2.get_start_coords().x - seg1.get_start_coords().x; //c1=a02
    system_extended_matrix[1][2] = seg2.get_start_coords().y - seg1.get_start_coords().y; //c2=a12
    system_extended_matrix[2][2] = seg2.get_start_coords().z - seg1.get_start_coords().z; //c3=a22

    //������� ��������� ����������� �������:
        //���� ��� =0 �� ������� ��������� ---�����
    for(int i=0; i<3; i++){
        for(int j=0; j<3;j++){
            if(system_extended_matrix[i][j] != 0 ) {  //����� ���������� ������ �� ���� � ����������� �����������, �� � �������� ������ ����������� �� �������
                    r_extend = 1;
                    break;
            }
        }
        if(r_extend==1) break;
    }
    if(r_extend == -1) { result.flag=2; return result; }

    //������� ������������ ����������� �������:
        //���� �� ����� ���� ������ ���� 3 --- ������� �����������, ������� �� ������������  -- �����
    extended_matrix_determinante = system_extended_matrix[0][0] * system_extended_matrix[1][1] * system_extended_matrix[2][2]
                                   + system_extended_matrix[1][0] *  system_extended_matrix[2][1] * system_extended_matrix[0][2]
                                   + system_extended_matrix[0][1] *  system_extended_matrix[1][2] * system_extended_matrix[2][0]
                                   - system_extended_matrix[0][2] *  system_extended_matrix[1][1] * system_extended_matrix[2][0]
                                   - system_extended_matrix[0][1] *  system_extended_matrix[1][0] * system_extended_matrix[2][2]
                                   - system_extended_matrix[1][2] *  system_extended_matrix[2][1] * system_extended_matrix[0][0];
    if(extended_matrix_determinante != 0) { result.flag=3; return result;  }

  //��������� ������ �������
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

    //���������� ������ ������� ������� � ���������� �� ����, � ������ ���������� ������ ��� ������� ������� ������� �������
    for(int i=0; i<3; i++){
       if( minors[i][0] != 0) { //����� ���������� ������ �� ���� � ����������� �����������, �� � �������� ������ ����������� �� �������
            r=2;
            delta_kramer = minors[i][0];
            delta_t_kramer = minors[i][1];
            delta_f_kramer = minors[i][2];
            break;
       }
       else r=1;
    }

    //���������� ������ ����������� ������� � ���������� �� ����
    for(int i=0; i<3; i++){
        for(int j=1;j<3;j++){
            if( minors[i][j] != 0) { r_extend=2; break; } //����� ���������� ������ �� ���� � ����������� �����������, �� � �������� ������ ����������� �� �������
            else r_extend=1;
        }
        if(r_extend == 2) break;
    }
    r_extend = max (r, r_extend);

    //���������� r � r_extend ���� �� ����� �� ������� �� ������������ -- �����
    if(r != r_extend) { result.flag=3; return result; }

    //���� ����� �� ������ 2 --- ������� ��������� --- �����
    if(r_extend <2 ) { result.flag=1; return result; }

    //���� ����� � ����� 2 �� ������ ������� ������� ������� ��������� ������������ ������ ������������
    //���� ����� ������ ����� -- ������ �������
    if(r_extend != 2) { result.flag=4; return result; }

    t= delta_t_kramer / delta_kramer;
    f= delta_f_kramer / delta_kramer;

    //�������� ��� seg1 ���������� ����� �� ���, ������������� ���������� t
    intersect_point_seg1=seg1.get_point_coords_defined_by_parameter(t);

    //�������� ��� seg2 ���������� ����� �� ���, ������������� ���������� f
    intersect_point_seg2=seg2.get_point_coords_defined_by_parameter(f);

    //���������� ���������� ���������� ����� ����������� �������� seg1 seg 2, ��� ������ ���������
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


// __******************************************** ������� ��� �������� ������� ********************************************
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


