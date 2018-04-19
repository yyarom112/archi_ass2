#include <stdio.h>
#include <stdlib.h>
extern void cmplx_add_s(double* a_real, double* a_img, double* b_real, double* b_img, double *res_real, double *res_img);
extern void cmplx_div_s(double* a_real, double* a_img, double* b_real, double* b_img, double *res_real, double *res_img);
extern void cmplx_sub_s(double* a_real, double* a_img, double* b_real, double* b_img, double *res_real, double *res_img);
extern void  cmplx_mult_s(double* a_real, double* a_img, double* b_real, double* b_img, double *res_real, double *res_img);
extern void eval_derivative_s(double *real_arr , double* img_arr,double* res_real_arr,double* res_img_arr,int len);
extern int* func_malloc();
extern void test_arr_float(double* arr,int len,double i);
void cmplx_add(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img);

void cmplx_sub(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img);

void cmplx_mult(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img);

void cmplx_div(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img);

void eval_poly(double *real_arr , double* img_arr, double *z_real, double* z_img, double *res_real, double *res_img, int pow);

void eval_derivative(double *real_arr , double* img_arr,double* res_real_arr,double* res_img_arr,int len);

double sqroot(double n);


double make_normal(double *a_real, double *a_img){
    double output,real=*a_real,img=*a_img;
    img*=img;
    real*=real;
    output=real+img;
    return sqroot(output);
}

//0=close enough 1=not close
int close_enough(double *a_real, double *a_img,double epsilon){
    if(make_normal(a_real, a_img)<epsilon)
        return 1;
    else
        return 0;
}

double sqroot(double n) {

    double temp=0, sqrt=0;
    sqrt=n/2;
    while(sqrt!=temp)
    {
        temp=sqrt;
        sqrt=(n/sqrt+sqrt)/2;
    }
    return sqrt;
}

void eval_poly(double *real_arr , double* img_arr, double *z_real, double* z_img, double *res_real, double *res_img, int pow){
    int i=1;
    *res_img=img_arr[0];
    *res_real=real_arr[0];
    double tmp_real=0,tmp_img=0;
    for(;i<=pow;i++){
        cmplx_mult(*res_real,*res_img,*z_real,*z_img,&tmp_real,&tmp_img);
        cmplx_add(tmp_real,tmp_img,real_arr[i],img_arr[i],res_real,res_img);
        tmp_real=0;
        tmp_img=0;
    }
}


void eval_derivative(double *real_arr , double* img_arr,double* res_real_arr,double* res_img_arr,int len){
    //x^2 + 2x+3
    int i=0,pow;
    for(;i<len;i++){
        pow=len-i;
        res_real_arr[i]=real_arr[i]*pow;
        res_img_arr[i]=img_arr[i]*pow;
    }
}

void cmplx_add(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img) {
    *res_real = a_real + b_real;
    *res_img =a_img+b_img;
}

void cmplx_mult(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img) {
    *res_real = a_real * b_real;
    *res_real+= -(a_img * b_img) ;
    *res_img=a_real * b_img;
    *res_img+=a_img * b_real;

}


void cmplx_sub(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img) {
    *res_real = a_real - b_real;
    *res_img=a_img-b_img;
}


void cmplx_div(double a_real, double a_img, double b_real, double b_img, double *res_real, double *res_img) {
    double divisor=(b_real*b_real)+(b_img*b_img);
    printf("divisor: %lf\n",divisor);
    *res_real = (a_real * b_real)+(a_img * b_img);
    *res_real=*res_real/divisor;
    *res_img=(a_img * b_real)-(a_real * b_img);
    *res_img=*res_img/divisor;

}

void newton_rashford_impl(double *real_arr , double* img_arr, double epsilon,int order, double* init_real, double* init_img){
    double* div_real_arr=malloc((order-1)* sizeof(double));
    double* div_img_arr=malloc((order-1)* sizeof(double));
    eval_derivative(real_arr,img_arr,div_real_arr,div_img_arr,order);
    //todo- on assembly we might assign the numbers on registers
    double f_real,f_img,fd_real,fd_img,c_real,c_img;

    eval_poly(real_arr,img_arr,init_real,init_img,&f_real,&f_img,order);
    eval_poly(div_real_arr,div_img_arr,init_real,init_img,&fd_real,&fd_img,order-1);
    cmplx_div(f_real,f_img,fd_real,fd_img,&c_real,&c_img);
    while(make_normal(&c_real,&c_img) >= epsilon){
        eval_poly(real_arr,img_arr,init_real,init_img,&f_real,&f_img,order);
        eval_poly(div_real_arr,div_img_arr,init_real,init_img,&fd_real,&fd_img,order-1);
        cmplx_div(f_real,f_img,fd_real,fd_img,&c_real,&c_img);
        //x(i+1)=x(i)-c
        *init_real-=c_real;
        *init_img-=c_img;
    }
    printf("the result:  (%lf , %lf)\n",*init_real,*init_img);

}

int main() {
    double a_real=1.3, a_img=0.2, b_real=5.2, b_img=0.1,res_real=0,res_img=0;
//    cmplx_div(a_real,a_img,b_real,b_img,&res_real,&res_img);
//    printf("%lf + %lfi\n",res_real,res_img);
//
    double arr_real[3]={5.3,2.0,3.0};
    double arr_img[3]={0.3,1.0,0.0};
    double res_arr_real[2]={0.0,0.0},res_arr_img[2]={0.0,0.0};
    printf("before dep: (%lf+%lf)x^2 +(%lf+%lf)x+(%lf+%lf)\n\n",arr_real[0],arr_img[0],arr_real[1],arr_img[1],arr_real[2],arr_img[2]);
    eval_derivative(arr_real,arr_img,res_arr_real,res_arr_img,2);
    printf("after dep: (%lf+%lf)x+(%lf+%lf)\n",res_arr_real[0],res_arr_img[0],res_arr_real[1],res_arr_img[1]);

    double res_arr_real_s[2]={0.0,0.0},res_arr_img_s[2]={0.0,0.0};
    eval_derivative_s(arr_real,arr_img,res_arr_real_s,res_arr_img_s,2);
    printf("after dep s: (%lf+%lf)x+(%lf+%lf)\n",res_arr_real_s[0],res_arr_img_s[0],res_arr_real_s[1],res_arr_img_s[1]);

//
//
//    double z_real=1;
//    double z_img=0;
//      double res_real=0;
//    double res_img=0;
//    printf("before dep: (%lf+%lf)x^2 +(%lf+%lf)x+(%lf+%lf)\n\n",arr_real[0],arr_img[0],arr_real[1],arr_img[1],arr_real[2],arr_img[2]);
//    eval_poly(arr_real,arr_img,&z_real,&z_img,&res_real,&res_img,2);
//    printf("the poly value for z= (%lf) is (%lf+%lfi)\n",z_real,res_real,res_img);
//
//
//    printf("sqroot %lf\n", sqroot(1.5));
//
//    printf("div %lf\n",make_normal(&res_real,&res_img));
//
//
//    printf("close_enoughs %d\n",close_enough(&res_real,&res_img,0.1));
//
//    double test1_real_arr[3]={1.0 , -4.0 , 4.0};
//    double test1_img_arr[3]={0.0 , 0.0 , 0.0};
//    double epsilon=0.00000000001;
//    int order=2;
//    double init_real=1.0, init_img=-1.0;
//    newton_rashford_impl(test1_real_arr ,test1_img_arr, epsilon, order,  &init_real,&init_img);
    printf("a: real=%lf, img=%lf\n",a_real,a_img);
    printf("b: real=%lf, img=%lf\n",b_real,b_img);

//    cmplx_add_s(&a_real,&a_img,&b_real,&b_img, &res_real, &res_img);
//
//    printf("add: real=%lf, img=%lf\n",res_real,res_img);
//    cmplx_sub_s(&a_real,&a_img,&b_real,&b_img, &res_real, &res_img);
//    printf("sub: real=%lf, img=%lf\n",res_real,res_img);
//    cmplx_mult_s(&a_real,&a_img,&b_real,&b_img, &res_real, &res_img);
//    printf("mult: real=%lf, img=%lf\n",res_real,res_img);



    cmplx_div(a_real,a_img,b_real,b_img, &res_real, &res_img);
    printf("div: real=%lf, img=%lf\n",res_real,res_img);

//    cmplx_div_s(&a_real,&a_img,&b_real,&b_img, &res_real, &res_img);
    printf("div_s: real=%lf, img=%lf\n",res_real,res_img);
    int *test=NULL;
    test= func_malloc();
    printf("div_s: real=%d,img:%d \n",test[0],test[1]);

//
//    double arr_real[3]={1.0,2.0,3.0};
//    double arr_img[3]={0.5,1.0,0.0};
//    double res_arr_real[2]={0.0,0.0},res_arr_img[2]={0.0,0.0};
//    printf("before dep: (%lf+%lf)x^2 +(%lf+%lf)x+(%lf+%lf)\n\n",arr_real[0],arr_img[0],arr_real[1],arr_img[1],arr_real[2],arr_img[2]);
//    eval_derivative(arr_real,arr_img,res_arr_real,res_arr_img,2);
//    printf("after dep: (%lf+%lf)x+(%lf+%lf)\n",res_arr_real[0],res_arr_img[0],res_arr_real[1],res_arr_img[1]);
//    double arr[2];
//    printf("before: arr[0]: %lf arr[1]:%lf\n",arr[0],arr[1]);
//
//    test_arr_float(arr,2,0.2);
//    printf("after arr[0]: %lf arr[1]:%lf\n",arr[0],arr[1]);

    return 0;
}
