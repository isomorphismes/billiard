#include "polygon.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static billiard_exact number(int32_t rational, int32_t radical, int32_t denominator) {
    billiard_exact result;
    assert(!billiard_exact_make(rational, radical, denominator, &result));
    return result;
}
static billiard_exact_point integral(int32_t x, int32_t y) {
    return (billiard_exact_point){number(x,0,1), number(y,0,1)};
}
static int same_number(billiard_exact first, billiard_exact second) {
    return first.rational == second.rational && first.radical == second.radical &&
           first.denominator == second.denominator;
}
static int same_ray(billiard_polygon_ray first, billiard_polygon_ray second) {
    return same_number(first.position.x, second.position.x) &&
           same_number(first.position.y, second.position.y) &&
           same_number(first.direction.x, second.direction.x) &&
           same_number(first.direction.y, second.direction.y) &&
           first.previous_side == second.previous_side &&
           first.steps_remaining == second.steps_remaining &&
           first.distance == second.distance && first.stopped == second.stopped &&
           first.field == second.field && first.launch_vertex == second.launch_vertex;
}
static void rational_equal(billiard_exact value, int64_t numerator, int64_t denominator) {
    assert(!value.radical);
    assert((int64_t)value.rational * denominator == numerator * value.denominator);
}

static billiard_polygon fixture(unsigned kind) {
    billiard_polygon polygon = {0};
    const int32_t coordinates[4][6][2] = {
        {{0,0},{1,0},{1,1},{0,1}},
        {{0,0},{5,0},{5,3},{0,3}},
        {{0,0},{2,0},{0,2}},
        {{0,0},{2,0},{2,1},{1,1},{1,2},{0,2}}
    };
    polygon.count = kind == 2 ? 3 : (kind == 3 ? 6 : 4);
    for (size_t vertex = 0; vertex < polygon.count; ++vertex) {
        polygon.vertex[vertex] = integral(coordinates[kind][vertex][0],
                                           coordinates[kind][vertex][1]);
        polygon.angle[vertex] = (billiard_angle){1,2};
    }
    if (kind == 2) polygon.angle[1] = polygon.angle[2] = (billiard_angle){1,4};
    if (kind == 3) polygon.angle[3] = (billiard_angle){3,2};
    return polygon;
}

/* Independent one-event oracle: integer line determinants at twice the
 * fixture coordinates. It does not call production field arithmetic,
 * polygon location, reflection, or collision selection. Small grid bounds
 * keep all signed 64-bit products below 10^8. */
typedef struct { size_t side, vertex; int64_t numerator, denominator; } oracle_hit;
static int64_t determinant(int64_t ax, int64_t ay, int64_t bx, int64_t by) {
    return ax*by-ay*bx;
}
static oracle_hit first_intersection(const billiard_polygon *polygon,
    int32_t start_x_twice, int32_t start_y_twice, int32_t dx, int32_t dy) {
    oracle_hit answer = {SIZE_MAX, SIZE_MAX, 0, 1};
    for (size_t side = 0; side < polygon->count; ++side) {
        size_t next = (side+1)%polygon->count;
        int64_t sx = 2*polygon->vertex[side].x.rational;
        int64_t sy = 2*polygon->vertex[side].y.rational;
        int64_t ex = 2*polygon->vertex[next].x.rational-sx;
        int64_t ey = 2*polygon->vertex[next].y.rational-sy;
        int64_t denominator = determinant(2*dx,2*dy,ex,ey);
        if (!denominator) continue;
        int64_t ray_numerator = determinant(sx-start_x_twice,sy-start_y_twice,ex,ey);
        int64_t side_numerator = determinant(sx-start_x_twice,sy-start_y_twice,2*dx,2*dy);
        if (denominator < 0) {
            denominator = -denominator; ray_numerator = -ray_numerator;
            side_numerator = -side_numerator;
        }
        if (ray_numerator <= 0 || side_numerator < 0 || side_numerator > denominator) continue;
        if (answer.side != SIZE_MAX &&
            ray_numerator*answer.denominator >= answer.numerator*denominator) continue;
        answer.side = side;
        answer.vertex = side_numerator == 0 ? side : (side_numerator == denominator ? next : SIZE_MAX);
        answer.numerator = ray_numerator; answer.denominator = denominator;
    }
    assert(answer.side != SIZE_MAX);
    return answer;
}

static void compare_grid(void) {
    unsigned compared = 0;
    for (unsigned kind = 0; kind < 4; ++kind) {
        billiard_polygon polygon = fixture(kind);
        assert(billiard_polygon_validate(&polygon) == BILLIARD_POLYGON_OK);
        int32_t xmax = kind == 0 ? 2 : (kind == 1 ? 10 : 4);
        int32_t ymax = kind == 0 ? 2 : (kind == 1 ? 6 : 4);
        for (int32_t x = 1; x < xmax; x += 2)
            for (int32_t y = 1; y < ymax; y += 2) {
                /* Closed-form membership, independent of winding-number code. */
                if ((kind == 2 && x+y >= 4) || (kind == 3 && x > 2 && y > 2)) continue;
                billiard_exact_point start = {number(x,0,2), number(y,0,2)};
                for (int32_t dx = -3; dx <= 3; ++dx)
                    for (int32_t dy = -3; dy <= 3; ++dy) {
                        if (!dx && !dy) continue;
                        oracle_hit expected = first_intersection(&polygon,x,y,dx,dy);
                        billiard_polygon_ray ray;
                        billiard_polygon_event event;
                        assert(billiard_polygon_begin(&polygon,start,integral(dx,dy),1,&ray) == BILLIARD_POLYGON_OK);
                        enum billiard_polygon_status wanted = expected.vertex == SIZE_MAX ?
                            BILLIARD_POLYGON_EDGE : BILLIARD_POLYGON_VERTEX_STOP;
                        assert(billiard_polygon_next(&polygon,&ray,&event) == wanted);
                        assert(event.side == expected.side && event.vertex == expected.vertex);
                        rational_equal(event.parameter,expected.numerator,expected.denominator);
                        rational_equal(event.end.x,x*expected.denominator+2*dx*expected.numerator,
                                       2*expected.denominator);
                        rational_equal(event.end.y,y*expected.denominator+2*dy*expected.numerator,
                                       2*expected.denominator);
                        if (expected.vertex == SIZE_MAX) {
                            size_t next = (expected.side+1)%polygon.count;
                            int64_t ex = polygon.vertex[next].x.rational-polygon.vertex[expected.side].x.rational;
                            int64_t ey = polygon.vertex[next].y.rational-polygon.vertex[expected.side].y.rational;
                            int64_t norm = ex*ex+ey*ey, projection = dx*ex+dy*ey;
                            rational_equal(event.outgoing.x,2*projection*ex-dx*norm,norm);
                            rational_equal(event.outgoing.y,2*projection*ey-dy*norm,norm);
                        }
                        double expected_length = hypot((double)dx,(double)dy)*
                            (double)expected.numerator/expected.denominator;
                        assert(fabs(event.length-expected_length) <= 1e-12);
                        billiard_polygon_ray saved = ray;
                        enum billiard_polygon_status stopped = expected.vertex == SIZE_MAX ?
                            BILLIARD_POLYGON_STEP_LIMIT : BILLIARD_POLYGON_STOPPED;
                        assert(billiard_polygon_next(&polygon,&ray,&event) == stopped);
                        assert(same_ray(ray,saved));
                        compared++;
                    }
            }
    }
    printf("independent integer collision/reflection oracle: %u cases PASS\n",compared);
}

static void equilateral(void) {
    billiard_polygon polygon = {0};
    polygon.count = 3;
    polygon.vertex[0] = integral(0,0); polygon.vertex[1] = integral(2,0);
    polygon.vertex[2] = (billiard_exact_point){number(1,0,1),number(0,1,1)};
    for (size_t vertex = 0; vertex < 3; ++vertex) polygon.angle[vertex] = (billiard_angle){1,3};
    assert(billiard_polygon_validate(&polygon) == BILLIARD_POLYGON_OK);
    billiard_exact_point start = {number(1,0,1),number(0,1,3)};
    billiard_polygon_ray ray;
    billiard_polygon_event event;
    assert(billiard_polygon_begin(&polygon,start,integral(1,0),4,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&polygon,&ray,&event) == BILLIARD_POLYGON_EDGE);
    assert(event.side == 1);
    assert(same_number(event.parameter,number(2,0,3)));
    assert(same_number(event.end.x,number(5,0,3)));
    assert(same_number(event.outgoing.x,number(-1,0,2)));
    assert(same_number(event.outgoing.y,number(0,-1,2)));
    assert(billiard_polygon_next(&polygon,&ray,&event) == BILLIARD_POLYGON_EDGE);
    assert(event.side == 0 && same_number(event.end.x,number(4,0,3)));
    assert(billiard_polygon_begin(&polygon,start,integral(0,1),4,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&polygon,&ray,&event) == BILLIARD_POLYGON_VERTEX_STOP);
    assert(event.vertex == 2 && same_number(event.end.y,number(0,1,1)));
    puts("equilateral exact sqrt(3) edge/vertex controls: PASS");
}

static billiard_polygon eighth_triangle_fixture(void) {
    billiard_polygon triangle = {0};
    triangle.field = BILLIARD_SQRT2_EIGHTHS;
    triangle.count = 3;
    triangle.vertex[0] = integral(0,0);
    triangle.vertex[1] = integral(1,0);
    triangle.vertex[2] = (billiard_exact_point){number(1,0,1),number(-1,1,1)};
    triangle.angle[0] = (billiard_angle){1,8};
    triangle.angle[1] = (billiard_angle){1,2};
    triangle.angle[2] = (billiard_angle){3,8};
    return triangle;
}

/* The diagonal is y=(sqrt(2)-1)x. For vertical incidence, its double
 * angle is pi/4, so the reflected vector is (sqrt(2)/2,-sqrt(2)/2).
 * This oracle uses closed-form coordinates independently of the kernel. */
static void eighth_angle_triangle(void) {
    billiard_polygon triangle = eighth_triangle_fixture();
    assert(billiard_polygon_validate(&triangle) == BILLIARD_POLYGON_OK);

    billiard_polygon_ray ray;
    billiard_polygon_event hit;
    for (int x = 8; x < 16; ++x) {
        billiard_exact_point start = {number(x,0,16),number(1,0,8)};
        assert(billiard_polygon_begin(&triangle,start,integral(0,1),2,&ray) == BILLIARD_POLYGON_OK);
        assert(billiard_polygon_next(&triangle,&ray,&hit) == BILLIARD_POLYGON_EDGE);
        assert(hit.side == 2 && hit.vertex == SIZE_MAX);
        assert(same_number(hit.parameter,number(-x-2,x,16)));
        assert(same_number(hit.end.x,start.x));
        assert(same_number(hit.end.y,number(-x,x,16)));
        assert(same_number(hit.outgoing.x,number(0,1,2)));
        assert(same_number(hit.outgoing.y,number(0,-1,2)));
        assert(fabs(hit.length-((sqrt(2.0)-1.0)*x/16.0-0.125)) < 1e-12);
        if (x == 8) {
            assert(billiard_polygon_next(&triangle,&ray,&hit) == BILLIARD_POLYGON_EDGE);
            assert(hit.side == 0 && same_number(hit.end.x,number(0,1,2)));
            assert(same_number(hit.parameter,number(2,-1,2)));
        }
    }
    billiard_exact_point start = {number(1,0,2),number(1,0,8)};
    assert(billiard_polygon_begin(&triangle,start,integral(-4,-1),2,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&triangle,&ray,&hit) == BILLIARD_POLYGON_VERTEX_STOP);
    assert(hit.vertex == 0 && same_number(hit.parameter,number(1,0,8)));
    assert(billiard_polygon_next(&triangle,&ray,&hit) == BILLIARD_POLYGON_STOPPED);
    /* 3-2*sqrt(2) is positive and below (sqrt(2)-1)/2. Using a sqrt(3)
     * sign predicate would place this point below the horizontal edge. */
    int location;
    billiard_exact_point near_edge = {number(1,0,2),number(3,-2,1)};
    assert(billiard_polygon_locate(&triangle,near_edge,&location) == BILLIARD_POLYGON_OK);
    assert(location == 1);

    /* Apply the similarity z -> (2+i)z. The first side is no longer on a
     * certified absolute heading; the relative angles are unchanged. */
    billiard_polygon rotated = triangle;
    rotated.vertex[1] = integral(2,1);
    rotated.vertex[2] = (billiard_exact_point){number(3,-1,1),number(-1,2,1)};
    assert(billiard_polygon_validate(&rotated) == BILLIARD_POLYGON_OK);
    billiard_polygon wrong = triangle;
    wrong.field = BILLIARD_SQRT3_TWELFTHS;
    assert(billiard_polygon_validate(&wrong) == BILLIARD_POLYGON_UNSUPPORTED_ANGLE);
    wrong = triangle;
    wrong.angle[0] = (billiard_angle){1,7};
    assert(billiard_polygon_validate(&wrong) == BILLIARD_POLYGON_UNSUPPORTED_ANGLE);
    wrong.field = (enum billiard_polygon_field)17;
    assert(billiard_polygon_validate(&wrong) == BILLIARD_POLYGON_INVALID);
    billiard_polygon square = fixture(0);
    assert(billiard_polygon_begin(&square,start,integral(0,1),2,&ray) == BILLIARD_POLYGON_OK);
    square.field = BILLIARD_SQRT2_EIGHTHS;
    assert(billiard_polygon_next(&square,&ray,&hit) == BILLIARD_POLYGON_INVALID);
    puts("exact pi/8 triangle intersections, reflections and vertex stop: PASS");
}

static void adverse(void) {
    billiard_polygon square = fixture(0), changed = square;
    billiard_exact_point normal;
    assert(billiard_polygon_side_normal(&square,0,&normal) == BILLIARD_POLYGON_OK);
    rational_equal(normal.x,0,1); rational_equal(normal.y,-1,1);
    assert(billiard_polygon_side_normal(&square,4,&normal) == BILLIARD_POLYGON_INVALID);
    billiard_polygon elbow = fixture(3);
    int location;
    assert(billiard_polygon_locate(&elbow,integral(1,1),&location) == BILLIARD_POLYGON_OK && location == 0);
    billiard_exact_point notch = {number(3,0,2),number(3,0,2)};
    assert(billiard_polygon_locate(&elbow,notch,&location) == BILLIARD_POLYGON_OK && location == -1);
    changed.vertex[1] = square.vertex[2]; changed.vertex[2] = square.vertex[1];
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);
    changed = square; changed.vertex[1] = changed.vertex[0];
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);
    changed = square; changed.vertex[1] = square.vertex[3]; changed.vertex[3] = square.vertex[1];
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);
    changed = square;
    for (size_t vertex = 0; vertex < 4; ++vertex) changed.vertex[vertex].y = number(0,0,1);
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);
    changed = square; changed.angle[0] = (billiard_angle){0,0};
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_UNSUPPORTED_ANGLE);
    changed = fixture(2); changed.vertex[1] = integral(3,0);
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_UNSUPPORTED_ANGLE);
    changed = square; changed.vertex[1].x = (billiard_exact){INT32_MIN,0,1};
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);
    changed = square; changed.count = BILLIARD_POLYGON_MAX_VERTICES+1;
    assert(billiard_polygon_validate(&changed) == BILLIARD_POLYGON_INVALID);

    billiard_polygon_ray ray;
    billiard_polygon_event event;
    billiard_exact_point center = {number(1,0,2),number(1,0,2)};
    assert(billiard_polygon_begin(&square,center,integral(1,1),4,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&square,&ray,&event) == BILLIARD_POLYGON_VERTEX_STOP);
    assert(event.vertex == 2); /* no reversal or continuation through a billiard vertex */
    billiard_exact_point almost_center = {number(1,0,2),number(500001,0,1000000)};
    assert(billiard_polygon_begin(&square,almost_center,integral(1,1),4,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&square,&ray,&event) == BILLIARD_POLYGON_EDGE);
    assert(event.side == 2 && event.vertex == SIZE_MAX);
    rational_equal(event.end.x,999999,1000000);
    billiard_exact_point almost_horizontal = {number(1,0,1),number(1,0,1000000)};
    assert(billiard_polygon_begin(&square,center,almost_horizontal,4,&ray) == BILLIARD_POLYGON_OK);
    assert(billiard_polygon_next(&square,&ray,&event) == BILLIARD_POLYGON_EDGE);
    assert(event.side == 1);
    rational_equal(event.end.y,1000001,2000000);
    assert(billiard_polygon_begin(&square,integral(0,0),integral(1,1),4,&ray) == BILLIARD_POLYGON_INVALID);
    assert(billiard_polygon_begin(&square,center,integral(0,0),4,&ray) == BILLIARD_POLYGON_INVALID);

    billiard_exact_point huge_direction = {number(INT32_MAX,0,1),number(INT32_MAX-1,0,1)};
    assert(billiard_polygon_begin(&square,center,huge_direction,4,&ray) == BILLIARD_POLYGON_OK);
    billiard_polygon_ray saved = ray;
    assert(billiard_polygon_next(&square,&ray,&event) == BILLIARD_POLYGON_ARITHMETIC_LIMIT);
    assert(same_ray(ray,saved));
    puts("malformed geometry, exact/near corner, nearly parallel, bounded arithmetic: PASS");
}

static void replay_and_adapter(void) {
    billiard_table table;
    assert(!billiard_default(&table,5,3));
    billiard_polygon rectangle;
    assert(billiard_polygon_rectangle(&table,&rectangle) == BILLIARD_POLYGON_OK);
    billiard_exact_point start = {number(5,0,2),number(1,0,1)};
    billiard_polygon_ray first, replay;
    assert(billiard_polygon_begin(&rectangle,start,integral(1,-1),80,&first) == BILLIARD_POLYGON_OK);
    replay = first;
    double distance = 0;
    for (unsigned step = 0; step < 80; ++step) {
        billiard_polygon_event event, repeated;
        enum billiard_polygon_status status = billiard_polygon_next(&rectangle,&first,&event);
        assert(status == BILLIARD_POLYGON_EDGE);
        assert(billiard_polygon_next(&rectangle,&replay,&repeated) == status);
        assert(event.side == repeated.side && same_number(event.parameter,repeated.parameter));
        assert(same_number(event.end.x,repeated.end.x) && same_number(event.end.y,repeated.end.y));
        distance += event.length;
        assert(first.distance == distance && replay.distance == distance);
        assert(fabs(billiard_fold(2.5+distance/sqrt(2.0),5)-billiard_exact_display(event.end.x)) < 1e-11);
        assert(fabs(billiard_fold(1-distance/sqrt(2.0),3)-billiard_exact_display(event.end.y)) < 1e-11);
    }
    billiard_polygon_event event;
    assert(billiard_polygon_next(&rectangle,&first,&event) == BILLIARD_POLYGON_STEP_LIMIT);
    table.width = NAN;
    assert(billiard_polygon_rectangle(&table,&rectangle) == BILLIARD_POLYGON_INVALID);
    table.width = 0x1p-40;
    assert(billiard_polygon_rectangle(&table,&rectangle) == BILLIARD_POLYGON_ARITHMETIC_LIMIT);
    puts("rectangle adapter, 80-event replay and independent mirror folding: PASS");
}

static void write_coefficients(billiard_exact value) {
    printf(" %" PRId32 " %" PRId32 " %" PRId32, value.rational,value.radical,value.denominator);
}

/* Emit exact results for the independently generated grid of interior rays.
 * Input consists of x*16, y*16 and the two integer direction components. */
static int trace_eighth_triangle(void) {
    billiard_polygon polygon = eighth_triangle_fixture();
    int x, y, dx, dy, count;
    while ((count = scanf("%d %d %d %d", &x,&y,&dx,&dy)) == 4) {
        billiard_polygon_ray ray;
        billiard_polygon_event hit;
        billiard_exact_point start = {number(x,0,16),number(y,0,16)};
        if (billiard_polygon_begin(&polygon,start,integral(dx,dy),1,&ray) != BILLIARD_POLYGON_OK)
            return 1;
        enum billiard_polygon_status status = billiard_polygon_next(&polygon,&ray,&hit);
        if (status != BILLIARD_POLYGON_EDGE && status != BILLIARD_POLYGON_VERTEX_STOP) return 1;
        printf("%u %zu", (unsigned)status, status == BILLIARD_POLYGON_EDGE ? hit.side : hit.vertex);
        write_coefficients(hit.parameter);
        write_coefficients(hit.end.x);
        write_coefficients(hit.end.y);
        write_coefficients(hit.outgoing.x);
        write_coefficients(hit.outgoing.y);
        putchar('\n');
    }
    return count == EOF && !ferror(stdin) && !ferror(stdout) ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1],"--trace-eighths")) return trace_eighth_triangle();
    if (argc != 1) return 1;
    compare_grid();
    equilateral();
    eighth_angle_triangle();
    adverse();
    replay_and_adapter();
    puts("BIL-POLY-01 bounded certified kernel: PASS");
    return 0;
}
