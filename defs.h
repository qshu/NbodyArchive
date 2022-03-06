typedef struct {
    short type;                 /* code for node type */
    real mass;                  /* total mass of node */
    vector pos;			/* position of node */
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)

/*
 * BODY: data structure used to represent particles.
 */

#define BODY 01                 /* type code for bodies */

typedef struct {
    short type;
    real mass;                  /* mass of body */
    vector pos;                 /* position of body */
    vector vel;                 /* velocity of body */
    vector acc;			/* acceleration of body */
    real phi;			/* potential at body */
} body, *bodyptr;

#define Body    body
#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)

/*
 * PHASEBODY: alternate definition introduced for I/O.
 */

typedef struct {
    short type;
    real mass;
    vector phase[2];            /* position, velocity of body */
    vector acc;
    real phi;
} phasebody;

#define Phase(x)  (((phasebody *) (x))->phase)

/*
 * CELL: structure used to represent internal nodes of tree.
 */

#define CELL 02                 /* type code for cells */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct {
    short type;
    real mass;                  /* total mass of cell */
    vector pos;                 /* cm. position of cell */
#ifdef QUADPOLE
    matrix quad;		/* quad. moment of cell */
#endif
    nodeptr subp[NSUB];         /* descendents of cell */
} cell, *cellptr;
