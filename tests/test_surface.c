#include "surface.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

/* Independent corner-cycle traversal: each corner has a neighbor across each
 * of its two incident sides. No production union/find or genus code is used. */
static void check_corner_cycles(const billiard_surface *surface) {
    unsigned char visited[BILLIARD_SURFACE_MAX_CORNERS] = {0};
    size_t queue[BILLIARD_SURFACE_MAX_CORNERS];
    size_t corners = surface->copy_count*surface->polygon.count, classes = 0;
    unsigned zero_sum = 0;
    for (size_t first = 0; first < corners; ++first) {
        if (visited[first]) continue;
        size_t begin = 0, end = 1;
        queue[0] = first; visited[first] = 1;
        unsigned angle_units = 0;
        while (begin < end) {
            size_t corner = queue[begin++];
            size_t copy = corner/surface->polygon.count, vertex = corner%surface->polygon.count;
            billiard_angle angle = surface->polygon.angle[vertex];
            assert((24u*angle.numerator)%angle.denominator == 0);
            angle_units += 24u*angle.numerator/angle.denominator;
            size_t incident[2] = {vertex,(vertex+surface->polygon.count-1)%surface->polygon.count};
            for (unsigned side = 0; side < 2; ++side) {
                size_t next = surface->seam[copy][incident[side]].target_copy*surface->polygon.count+vertex;
                if (!visited[next]) { visited[next] = 1; queue[end++] = next; }
            }
        }
        assert(angle_units%48 == 0);
        size_t label = surface->vertex_class[first/surface->polygon.count][first%surface->polygon.count];
        assert(surface->cone[label].angle_twelfths == angle_units/2);
        zero_sum += angle_units/48-1;
        ++classes;
    }
    assert(classes == surface->vertex_class_count);
    int euler = (int)classes-(int)(corners/2)+(int)surface->copy_count;
    assert(euler == surface->euler_characteristic);
    assert(2-(int)(2*surface->genus) == euler);
    assert(zero_sum == 2*surface->genus-2);
}

static void compare_files(FILE *first, FILE *second) {
    rewind(first); rewind(second);
    for (;;) {
        int left = fgetc(first), right = fgetc(second);
        assert(left == right);
        if (left == EOF) break;
    }
}

/* The dihedral group acts as z -> exp(i*k*pi/4) times z or its
 * conjugate. Read k from the image of (1,0), without using production
 * field operations or the generated reflection matrices. */
static unsigned octagonal_heading(const billiard_surface_copy *copy) {
    const exact_number cosine[8] = {
        {1,0,1},{0,1,2},{0,0,1},{0,-1,2},
        {-1,0,1},{0,-1,2},{0,0,1},{0,1,2}
    };
    for (unsigned k = 0; k < 8; ++k) {
        if (billiard_exact_equal(copy->placement.entry[0][0],cosine[k]) &&
            billiard_exact_equal(copy->placement.entry[1][0],cosine[(k+6)%8]))
            return k;
    }
    assert(0 && "unexpected image of horizontal unit vector");
    return 0;
}

static void eighth_triangle_cover(void) {
    billiard_polygon triangle = {0};
    triangle.field = BILLIARD_SQRT2_EIGHTHS;
    triangle.count = 3;
    triangle.vertex[0] = (exact_point){exact_integer(0),exact_integer(0)};
    triangle.vertex[1] = (exact_point){exact_integer(1),exact_integer(0)};
    triangle.vertex[2] = (exact_point){exact_integer(1),(exact_number){-1,1,1}};
    triangle.angle[0] = (billiard_angle){1,8};
    triangle.angle[1] = (billiard_angle){1,2};
    triangle.angle[2] = (billiard_angle){3,8};

    billiard_surface surface;
    assert(billiard_surface_unfold(&triangle,&surface) == BILLIARD_SURFACE_OK);
    assert(billiard_surface_validate(&surface) == BILLIARD_SURFACE_OK);
    assert(billiard_surface_validate_twelfths(&surface) == BILLIARD_SURFACE_UNSUPPORTED_FIELD);
    assert(surface.copy_count == 16 && surface.edge_count == 24 && surface.vertex_class_count == 6);
    assert(surface.components == 1 && surface.euler_characteristic == -2 && surface.genus == 2);
    assert(billiard_exact_equal(surface.area,(exact_number){-8,8,1}));
    check_corner_cycles(&surface);

    int seen[2][8] = {{0}};
    const int double_edge_angle[] = {0,4,1};
    exact_work work = {.field = triangle.field};
    for (size_t copy = 0; copy < surface.copy_count; ++copy) {
        unsigned heading = octagonal_heading(&surface.copy[copy]);
        int orientation = surface.copy[copy].orientation;
        assert(orientation == 1 || orientation == -1);
        assert(!seen[orientation<0][heading]++);
        for (size_t side = 0; side < 3; ++side) {
            size_t target = surface.seam[copy][side].target_copy;
            assert(target < 16 && surface.copy[target].orientation == -orientation);
            unsigned expected = (unsigned)(8+(int)heading+orientation*double_edge_angle[side])%8;
            assert(octagonal_heading(&surface.copy[target]) == expected);
            exact_point vertex = billiard_surface_vertex(&work,&surface,copy,side);
            exact_point folded = billiard_surface_fold(&work,&surface,copy,vertex);
            assert(!work.failed && exact_point_equal(folded,triangle.vertex[side]));
        }
    }
    size_t vertex_counts[3] = {0};
    for (size_t number = 0; number < surface.vertex_class_count; ++number) {
        size_t source = surface.cone[number].representative % 3;
        ++vertex_counts[source];
        assert(surface.cone[number].zero_order == (source == 2 ? 2u : 0u));
    }
    assert(vertex_counts[0] == 1 && vertex_counts[1] == 4 && vertex_counts[2] == 1);

    FILE *written = tmpfile(); assert(written && !billiard_surface_write(written,&surface));
    rewind(written);
    char line[32];
    assert(fgets(line,sizeof line,written) && strcmp(line,"radicand\t2\n") == 0);
    fclose(written);
    surface.area.radical = -8;
    assert(billiard_surface_validate(&surface) == BILLIARD_SURFACE_INVALID);
    puts("pi/8 triangle: 16 dihedral charts, 48 seams, genus two, one order-two zero PASS");
}

int main(void) {
    eighth_triangle_cover();
    const unsigned copies[] = {4,4,8,6,4}, vertices[] = {4,4,4,3,6};
    for (unsigned fixture = 0; fixture < BILLIARD_FIXTURE_COUNT; ++fixture) {
        billiard_polygon polygon;
        billiard_surface surface;
        assert(billiard_polygon_fixture((enum billiard_fixture)fixture,&polygon) == BILLIARD_POLYGON_OK);
        assert(billiard_surface_unfold(&polygon,&surface) == BILLIARD_SURFACE_OK);
        assert(billiard_surface_validate(&surface) == BILLIARD_SURFACE_OK);
        assert(surface.copy_count == copies[fixture]);
        assert(surface.vertex_class_count == vertices[fixture]);
        assert(surface.genus == (fixture == BILLIARD_L_POLYGON ? 2u : 1u));
        if (fixture == BILLIARD_L_POLYGON) assert(billiard_exact_equal(surface.area,exact_integer(12)));
        if (fixture == BILLIARD_SQUARE) assert(billiard_exact_equal(surface.area,exact_integer(4)));
        if (fixture == BILLIARD_RECTANGLE) assert(billiard_exact_equal(surface.area,exact_integer(60)));
        if (fixture == BILLIARD_RIGHT_TRIANGLE) assert(billiard_exact_equal(surface.area,exact_integer(16)));
        if (fixture == BILLIARD_EQUILATERAL) assert(billiard_exact_equal(surface.area,(exact_number){0,6,1}));
        check_corner_cycles(&surface);
        exact_work work = {0};
        for (size_t copy = 0; copy < surface.copy_count; ++copy)
            for (size_t vertex = 0; vertex < polygon.count; ++vertex) {
                exact_point placed = billiard_surface_vertex(&work,&surface,copy,vertex);
                exact_point folded = billiard_surface_fold(&work,&surface,copy,placed);
                assert(exact_point_equal(folded,polygon.vertex[vertex]));
            }
        assert(!work.failed);
        FILE *first = tmpfile(), *second = tmpfile(); assert(first && second);
        assert(!billiard_surface_write(first,&surface) && !billiard_surface_write(second,&surface));
        compare_files(first,second); fclose(first); fclose(second);

        billiard_surface changed = surface;
        changed.seam[0][0].target_copy = SIZE_MAX;
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        changed = surface; changed.seam[0][0].translation.x = exact_integer(17);
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        changed = surface; changed.copy[1].orientation *= -1;
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        changed = surface; changed.copy_count--;
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        changed = surface; changed.cone[0].angle_twelfths++;
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        changed = surface; changed.vertex_class[0][0] = SIZE_MAX;
        assert(billiard_surface_validate(&changed) == BILLIARD_SURFACE_INVALID);
        printf("unfold fixture %u: copies=%zu V=%zu E=%u F=%zu genus=%u independent cycles PASS\n",
            fixture,surface.copy_count,surface.vertex_class_count,surface.edge_count,surface.copy_count,surface.genus);
    }
    billiard_polygon polygon;
    billiard_surface surface;
    assert(billiard_polygon_fixture(BILLIARD_L_POLYGON,&polygon) == BILLIARD_POLYGON_OK);
    polygon.field = BILLIARD_SQRT2_EIGHTHS;
    assert(billiard_surface_unfold(&polygon,&surface) == BILLIARD_SURFACE_OK);
    assert(surface.genus == 2);
    polygon.field = BILLIARD_SQRT3_TWELFTHS;
    polygon.angle[3] = (billiard_angle){0,0};
    assert(billiard_surface_unfold(&polygon,&surface) == BILLIARD_SURFACE_INVALID);
    puts("BIL-UNFOLD-02 exact covers, topology and malformed-gluing controls: PASS");
    return 0;
}
