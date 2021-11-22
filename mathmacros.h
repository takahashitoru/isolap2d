#ifndef MATHMACROS_H
#define MATHMACROS_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef PI
#define PI (M_PI)
#else
#error PI is already defined.
#endif

#ifndef PI2I
#define PI2I 0.15915494309189533576 // 1/(2*PI)
#endif

#ifndef GAM
#define GAM 0.57721566490153286060
#else
#error GAM is already defined.
#endif

#ifndef MAX
#define MAX(a, b) ( (a) >= (b) ? (a) : (b) )
#else
#error MAX is already defined.
#endif

#ifndef MIN
#define MIN(a, b) ( (a) <= (b) ? (a) : (b) )
#else
#error MIN is already defined.
#endif

#ifndef SQUARE
#define SQUARE(x) ( (x) * (x) )
#else
#error SQUARE is already defined.
#endif

#ifndef CUBE
#define CUBE(x) ( (x) * (x) * (x) )
#else
#error CUBE is already defined.
#endif

#ifndef RAND
#define RAND(xmin, xmax) ( (double)(xmin) + (double)((xmax) - (xmin)) * rand() / (double)RAND_MAX )
#else
#error RAND is already defined.
#endif

#endif /* MATHMACROS_H */
