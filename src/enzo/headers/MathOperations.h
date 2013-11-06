//
// Mathematical Operations
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __MATH_OPERATIONS_H__
#define __MATH_OPERATIONS_H__

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

typedef float (*MathFunction)(float, float);

inline float MinVal(float a, float b) {
    if (a < b) return a;
    return b;
}

inline float MaxVal(float a, float b) {
    if (a > b) return a;
    return b;
}

inline float CopyVal(float a, float b) {
    return b;
}

inline float AddVal(float a, float b) {
    return (a + b);
}

inline float SubVal(float a, float b) {
    return (a - b);
}

inline float MultVal(float a, float b) {
    return a*b;
}

inline float DivVal(float a, float b) {
    return a/b;
}

#endif
