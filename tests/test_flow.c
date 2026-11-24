#include "flow.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static exact_number fraction(int32_t numerator, int32_t denominator) {
    exact_number result; assert(!billiard_exact_make(numerator,0,denominator,&result)); return result;
}
static exact_point integral(int32_t horizontal, int32_t vertical) {
    return (exact_point){exact_integer(horizontal),exact_integer(vertical)};
}
static billiard_surface fixture(enum billiard_fixture kind) {
    billiard_polygon polygon; billiard_surface surface;
    assert(billiard_polygon_fixture(kind,&polygon) == BILLIARD_POLYGON_OK);
    assert(billiard_surface_unfold(&polygon,&surface) == BILLIARD_SURFACE_OK);
    return surface;
}

static enum billiard_flow_status run(const billiard_surface *surface, billiard_flow *flow,
    double width, double height, unsigned *regular_vertices) {
    enum billiard_flow_status status;
    for (unsigned step = 0; step < 300; ++step) {
        billiard_flow_event event;
        status = billiard_flow_next(surface,flow,&event);
        if (status != BILLIARD_FLOW_EDGE && status != BILLIARD_FLOW_REGULAR_VERTEX &&
            status != BILLIARD_FLOW_PERIODIC && status != BILLIARD_FLOW_SINGULAR &&
            status != BILLIARD_FLOW_DISTANCE_LIMIT) return status;
        if (status == BILLIARD_FLOW_REGULAR_VERTEX) ++*regular_vertices;
        assert(exact_point_equal(event.global_direction,flow->direction));
        if (width > 0) {
            double parameter = billiard_exact_display(flow->parameter);
            double expected_x = billiard_fold(billiard_exact_display(flow->initial_position.x)+
                parameter*billiard_exact_display(flow->direction.x),width);
            double expected_y = billiard_fold(billiard_exact_display(flow->initial_position.y)+
                parameter*billiard_exact_display(flow->direction.y),height);
            assert(fabs(expected_x-billiard_exact_display(event.folded_end.x)) < 1e-10);
            assert(fabs(expected_y-billiard_exact_display(event.folded_end.y)) < 1e-10);
        }
        if (status != BILLIARD_FLOW_EDGE && status != BILLIARD_FLOW_REGULAR_VERTEX) return status;
    }
    return BILLIARD_FLOW_CROSSING_LIMIT;
}

static void periodic_control(enum billiard_fixture kind, exact_point position,
    exact_point direction, int32_t expected_period, double width, double height, int expect_vertex) {
    billiard_surface surface = fixture(kind);
    billiard_flow flow, replay;
    assert(billiard_flow_begin(&surface,0,position,direction,exact_integer(40),200,&flow) == BILLIARD_FLOW_OK);
    replay = flow;
    unsigned regular = 0, repeated_regular = 0;
    assert(run(&surface,&flow,width,height,&regular) == BILLIARD_FLOW_PERIODIC);
    assert(run(&surface,&replay,width,height,&repeated_regular) == BILLIARD_FLOW_PERIODIC);
    assert(billiard_exact_equal(flow.parameter,exact_integer(expected_period)));
    assert(billiard_exact_equal(flow.parameter,replay.parameter));
    assert(flow.copy == replay.copy && flow.crossings == replay.crossings && regular == repeated_regular);
    if (expect_vertex) assert(regular > 0);
    printf("flow fixture %u: exact first-return parameter %d, regular vertices %u PASS\n",kind,expected_period,regular);
}

static billiard_surface eighth_triangle(void) {
    billiard_polygon triangle = {0};
    triangle.field = BILLIARD_SQRT2_EIGHTHS;
    triangle.count = 3;
    triangle.vertex[0] = integral(0,0);
    triangle.vertex[1] = integral(1,0);
    triangle.vertex[2] = (exact_point){exact_integer(1),(exact_number){-1,1,1}};
    triangle.angle[0] = (billiard_angle){1,8};
    triangle.angle[1] = (billiard_angle){1,2};
    triangle.angle[2] = (billiard_angle){3,8};
    billiard_surface surface;
    assert(billiard_surface_unfold(&triangle,&surface) == BILLIARD_SURFACE_OK);
    return surface;
}

static void eighth_flow_control(void) {
    billiard_surface surface = eighth_triangle();
    exact_point start = {fraction(1,2),fraction(1,8)};
    billiard_flow flow;
    billiard_flow_event event;
    assert(billiard_flow_begin(&surface,0,start,integral(0,1),exact_integer(4),8,&flow) == BILLIARD_FLOW_OK);

    /* Intersect the folded lines with y=(sqrt(2)-1)x, y=0, x=1 and
     * y=(sqrt(2)-1)x again. The developed line stays vertical at x=1/2. */
    const struct {
        size_t side;
        exact_number time, folded_x, folded_y, developed_y;
    } expected[] = {
        {2, {-5,4,8}, {1,0,2}, {-1,1,2}, {-1,1,2}},
        {0, {2,-1,2}, {0,1,2}, {0,0,1},  {1,0,2}},
        {1, {-1,1,1}, {1,0,1},  {2,-1,2}, {-1,2,2}},
        {2, {-4,3,2}, {-1,2,2}, {5,-3,2}, {-5,5,2}}
    };
    for (size_t step = 0; step < sizeof expected/sizeof *expected; ++step) {
        assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
        assert(event.field == BILLIARD_SQRT2_EIGHTHS && event.side == expected[step].side);
        assert(event.to_copy != event.from_copy && flow.crossings == step+1);
        assert(billiard_exact_equal(event.parameter,expected[step].time));
        assert(billiard_exact_equal(event.folded_end.x,expected[step].folded_x));
        assert(billiard_exact_equal(event.folded_end.y,expected[step].folded_y));
        assert(billiard_exact_equal(event.developed_end.x,start.x));
        assert(billiard_exact_equal(event.developed_end.y,expected[step].developed_y));
        assert(exact_point_equal(event.global_direction,integral(0,1)));
        exact_number time = expected[step].time;
        double length = (time.rational+sqrt(2.0)*time.radical)/time.denominator;
        assert(fabs(event.length-length) < 1e-13);
    }
    FILE *output = tmpfile();
    assert(output && !billiard_flow_write_event(output,&event));
    rewind(output);
    char line[512];
    assert(fgets(line,sizeof line,output) && strstr(line,"\tradicand=2\n"));
    fclose(output);

    /* Launch the same folded ray from the reflected chart. Its matrix is
     * sqrt(2)/2 * [[1,1],[1,-1]], so both the global direction and the
     * developed start contain radicals before the first flight. */
    size_t reflected_chart = surface.seam[0][2].target_copy;
    exact_point lift_direction = {(exact_number){0,1,2},(exact_number){0,-1,2}};
    assert(billiard_flow_begin(&surface,reflected_chart,start,lift_direction,exact_integer(2),4,&flow)
           == BILLIARD_FLOW_OK);
    assert(exact_point_equal(flow.folded_ray.direction,integral(0,1)));
    exact_point lift_start = {(exact_number){0,5,16},(exact_number){0,3,16}};
    assert(exact_point_equal(flow.initial_chart_position,lift_start));
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
    assert(event.from_copy == reflected_chart && event.to_copy == 0 && event.side == 2);
    assert(exact_point_equal(event.developed_end,event.folded_end));
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
    assert(event.side == 0 && billiard_exact_equal(event.developed_end.x,(exact_number){0,1,2}));
    assert(billiard_exact_equal(event.developed_end.y,exact_integer(0)));

    /* Aim directly at the 3*pi/8 corner. Its 16 lifts meet at a 6*pi
     * cone, so continuing the same vector past this event is undefined. */
    exact_point to_cone = {fraction(1,2),(exact_number){-9,8,8}};
    assert(billiard_flow_begin(&surface,0,start,to_cone,exact_integer(2),8,&flow) == BILLIARD_FLOW_OK);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_SINGULAR);
    assert(event.vertex_class == surface.vertex_class[0][2] && surface.cone[event.vertex_class].zero_order == 2);
    assert(billiard_exact_equal(event.parameter,exact_integer(1)));
    assert(exact_point_equal(event.developed_end,surface.polygon.vertex[2]));
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_STOPPED);

    exact_point to_right_angle = {fraction(1,2),fraction(-1,8)};
    assert(billiard_flow_begin(&surface,0,start,to_right_angle,exact_integer(2),8,&flow) == BILLIARD_FLOW_OK);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_REGULAR_VERTEX);
    assert(surface.cone[event.vertex_class].zero_order == 0 && event.to_copy != event.from_copy);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_DISTANCE_LIMIT);
    assert(exact_point_equal(event.folded_end,start));
    assert(exact_point_equal(event.developed_end,(exact_point){fraction(3,2),fraction(-1,8)}));
    puts("pi/8 genus-two cover: four analytic seam events, regular continuation, singular stop PASS");
}

static void radical_width_return(void) {
    billiard_polygon rectangle;
    billiard_surface surface;
    assert(billiard_polygon_fixture(BILLIARD_SQUARE,&rectangle) == BILLIARD_POLYGON_OK);
    rectangle.field = BILLIARD_SQRT2_EIGHTHS;
    rectangle.vertex[1].x = rectangle.vertex[2].x = (exact_number){0,1,1};
    assert(billiard_surface_unfold(&rectangle,&surface) == BILLIARD_SURFACE_OK);
    exact_point start = {fraction(1,2),fraction(1,2)};
    billiard_flow flow;
    billiard_flow_event event;
    assert(billiard_flow_begin(&surface,0,start,integral(1,0),exact_integer(3),4,&flow) == BILLIARD_FLOW_OK);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_PERIODIC);
    assert(flow.copy == 0 && flow.crossings == 2 && exact_point_equal(event.folded_end,start));
    assert(billiard_exact_equal(flow.parameter,(exact_number){0,2,1}));
    puts("sqrt(2)-wide rectangle: exact first return at 2*sqrt(2) PASS");
}

int main(void) {
    eighth_flow_control();
    radical_width_return();
    exact_point square_start = {fraction(1,2),fraction(1,3)};
    periodic_control(BILLIARD_SQUARE,square_start,integral(1,0),2,1,1,0);
    periodic_control(BILLIARD_SQUARE,square_start,integral(2,1),2,1,1,0);
    exact_point center = {fraction(1,2),fraction(1,2)};
    periodic_control(BILLIARD_SQUARE,center,integral(1,1),2,1,1,1);
    periodic_control(BILLIARD_RECTANGLE,(exact_point){fraction(5,2),exact_integer(1)},integral(5,3),2,5,3,0);
    periodic_control(BILLIARD_RIGHT_TRIANGLE,center,integral(1,0),4,0,0,0);
    periodic_control(BILLIARD_EQUILATERAL,(exact_point){exact_integer(1),(exact_number){0,1,3}},integral(1,0),6,0,0,0);
    periodic_control(BILLIARD_L_POLYGON,center,integral(1,0),4,0,0,0);
    periodic_control(BILLIARD_L_POLYGON,(exact_point){fraction(1,2),fraction(3,2)},integral(1,0),2,0,0,0);

    billiard_surface surface = fixture(BILLIARD_L_POLYGON);
    billiard_flow flow;
    exact_point above_notch = {fraction(1,2),fraction(3,2)};
    assert(billiard_flow_begin(&surface,0,above_notch,integral(1,-1),exact_integer(10),40,&flow) == BILLIARD_FLOW_OK);
    billiard_flow_event event;
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_SINGULAR);
    assert(billiard_exact_equal(flow.parameter,fraction(1,2)));
    assert(surface.cone[event.vertex_class].zero_order == 2);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_STOPPED);

    surface = fixture(BILLIARD_SQUARE);
    exact_point irrational = {exact_integer(1),(exact_number){0,1,1}};
    assert(billiard_flow_begin(&surface,0,square_start,irrational,exact_integer(10),80,&flow) == BILLIARD_FLOW_OK);
    unsigned regular = 0;
    enum billiard_flow_status status = run(&surface,&flow,1,1,&regular);
    assert(status == BILLIARD_FLOW_DISTANCE_LIMIT || status == BILLIARD_FLOW_ARITHMETIC_LIMIT || status == BILLIARD_FLOW_CROSSING_LIMIT);
    exact_point almost_periodic = {exact_integer(1),fraction(1001,1000)};
    assert(billiard_flow_begin(&surface,0,square_start,almost_periodic,exact_integer(10),80,&flow) == BILLIARD_FLOW_OK);
    status = run(&surface,&flow,1,1,&regular);
    assert(status == BILLIARD_FLOW_DISTANCE_LIMIT || status == BILLIARD_FLOW_ARITHMETIC_LIMIT || status == BILLIARD_FLOW_CROSSING_LIMIT);
    assert(billiard_flow_begin(&surface,0,square_start,integral(1,0),exact_integer(10),1,&flow) == BILLIARD_FLOW_OK);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_EDGE);
    exact_number saved_parameter = flow.parameter;
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_CROSSING_LIMIT);
    assert(billiard_exact_equal(flow.parameter,saved_parameter));
    assert(billiard_flow_begin(&surface,0,square_start,integral(1,0),fraction(1,4),10,&flow) == BILLIARD_FLOW_OK);
    assert(billiard_flow_next(&surface,&flow,&event) == BILLIARD_FLOW_DISTANCE_LIMIT);
    assert(event.side == SIZE_MAX && billiard_exact_equal(event.folded_end.x,fraction(3,4)));
    assert(billiard_flow_begin(&surface,0,square_start,integral(0,0),exact_integer(10),10,&flow) == BILLIARD_FLOW_INVALID);
    assert(billiard_flow_begin_vertex(&surface,0,0,integral(1,1),exact_integer(10),40,&flow) == BILLIARD_FLOW_OK);
    regular = 0;
    assert(run(&surface,&flow,0,0,&regular) == BILLIARD_FLOW_PERIODIC);
    assert(billiard_exact_equal(flow.parameter,exact_integer(2)));
    assert(billiard_flow_begin_vertex(&surface,0,0,integral(-1,-1),exact_integer(10),40,&flow) == BILLIARD_FLOW_INVALID);
    assert(billiard_flow_begin(&surface,0,square_start,integral(2,1),exact_integer(10),40,&flow) == BILLIARD_FLOW_OK);
    billiard_flow replay = flow;
    FILE *first = tmpfile(), *second = tmpfile(); assert(first && second);
    for (unsigned step = 0; step < 40; ++step) {
        billiard_flow_event left, right;
        enum billiard_flow_status left_status = billiard_flow_next(&surface,&flow,&left);
        enum billiard_flow_status right_status = billiard_flow_next(&surface,&replay,&right);
        assert(left_status == right_status);
        assert(!billiard_flow_write_event(first,&left) && !billiard_flow_write_event(second,&right));
        if (left_status == BILLIARD_FLOW_PERIODIC) break;
        assert(left_status == BILLIARD_FLOW_EDGE || left_status == BILLIARD_FLOW_REGULAR_VERTEX);
    }
    rewind(first); rewind(second);
    for (;;) { int left = fgetc(first), right = fgetc(second); assert(left == right); if (left == EOF) break; }
    fclose(first); fclose(second);
    surface.seam[0][0].target_copy = SIZE_MAX;
    assert(billiard_flow_begin(&surface,0,square_start,integral(1,0),exact_integer(10),10,&flow) == BILLIARD_FLOW_INVALID);
    puts("BIL-FLOW-03 straight seams, independent mirror folding, regular vertices, singular stops and UNKNOWN limits: PASS");
    return 0;
}
