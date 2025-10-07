#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>

typedef union{
    int32_t i;
    float f;
    struct{
        uint32_t mantissa : 23;
        uint32_t expoent: 8;
        uint32_t sign: 1;
    } parts;
} Float_t;

void printFloat_t(Float_t num){
    printf("f:%1.9e, ix:0x%08X, s:%d, e:%d, mx:0x%06X\n\n",
    num.f, num.i, num.parts.sign, num.parts.expoent, num.parts.mantissa);
}

int main(){
    Float_t num;

    printf("Epsilon: %1.15f\n\n", FLT_EPSILON);

    num.f = 0.0f;
    printFloat_t(num);


    num.f = 1.40129846e-45f;
    printFloat_t(num);


    num.f = 1.17549435e-38f;
    printFloat_t(num);


    num.f = 0.2f;
    printFloat_t(num);


    num.f = 1.0f;
    printFloat_t(num);


    num.f = 1.5f;
    printFloat_t(num); 


    num.f = 1.75f;
    printFloat_t(num);


    num.f = 1.99999988f;
    printFloat_t(num);


    num.f = 2.0f;
    printFloat_t(num);


    num.f = 16777215.0f;
    printFloat_t(num);


    num.f = 3.40282347e+38f;
    printFloat_t(num);


    num.f = INFINITY;
    printFloat_t(num);

    return 0;
}
