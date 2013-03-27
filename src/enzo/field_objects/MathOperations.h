//
// Mathematical Operations
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __MATH_OPERATIONS_H__
#define __MATH_OPERATIONS_H__

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))

typedef double (*MathFunction)(double, double);

inline double MinVal(double a, double b) {
    if (a < b) return a;
    return b;
}

inline double MaxVal(double a, double b) {
    if (a > b) return a;
    return b;
}

inline double CopyVal(double a, double b) {
    return b;
}

inline double AddVal(double a, double b) {
    return (a + b);
}

inline double SubVal(double a, double b) {
    return (a - b);
}

inline double MultVal(double a, double b) {
    return a*b;
}

inline double DivVal(double a, double b) {
    return a/b;
}

#endif
