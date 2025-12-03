//as mentioned in the paper better to use to adjust it to different dimensions without code changes
#define DIM 2
typedef double vect_t[DIM];
struct Body
{
    double mass;
    vect_t s; // position
    vect_t v; // velocity
    vect_t a; // acceleration
    vect_t f; // force
};
