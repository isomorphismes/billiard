# Certified polygon collisions — BIL-POLY-01

The kernel follows a point in a polygon until its next boundary intersection.
It reflects the direction at a regular edge and stops at a vertex. Every
intersection, ordering comparison, angle check and reflection uses exact
arithmetic. A close pair of events is never merged using a tolerance.

This implementation is bounded. Coordinates have the form (a + b√r)/d, with r = 2 or 3,
32-bit integer coefficients and a positive denominator. Checked operations
return ARITHMETIC_LIMIT if an intermediate value cannot be represented. No
event or ray state is committed on failure. Conservative limits can occur
even when the eventual reduced answer would fit. This is not arbitrary
precision arithmetic or support for every rational-angle polygon.

The polygon selects its coordinate field and certified heading family. A
zero-initialized polygon uses Q(√3) with headings in multiples of π/12.
Set `field = BILLIARD_SQRT2_EIGHTHS` for Q(√2) and multiples of π/8. The
headings are relative to the first side in both cases. The eighth-angle
directions follow from tan(π/8) = √2 − 1, swapping coordinates and quarter
turns. Exact determinants check parallelism and dot products choose the
orientation. Validation verifies every claimed interior angle against the
vertices. Missing or incorrect claims and unsupported directions fail.

All vertices, rays and events of one polygon use its field. Use
`billiard_polygon_display(polygon.field, number)` for drawing radical
coordinates. The standalone exact arithmetic and translation-cover API
retain their Q(√3) contract. Unfolding an eighth-angle polygon returns
`BILLIARD_SURFACE_UNSUPPORTED_FIELD`.

The π/8 regression triangle has vertices (0,0), (1,0), (1,√2−1).
Vertical incidence at x hits the diagonal at y=(√2−1)x and leaves in
direction (√2/2, −√2/2). Eight starting positions check this closed form.

Geometry is counterclockwise in its coordinate plane. Validation rejects
zero edges, backtracking, nonadjacent crossings or touching, nonpositive
area and incorrect angle claims. Collinear forward edges are allowed with
an interior angle π. Point location uses exact winding and boundary tests.
The kernel receives no pixels, roughness, window, timing or touch state.

The unit square, nonsquare rectangle, right-isosceles triangle, equilateral
triangle with exact √3 height, and three-square right-angled L domain are
controls. A separate integer-determinant oracle checks every first event
and reflected direction on a bounded grid in the rational fixtures. The
equilateral controls use hand-derived exact intersections and reflected
vectors. Repeated rectangle events are compared with mirror folding.
Negative cases include self-crossing, clockwise, zero-area and untagged
input, an unsupported triangle, exact versus nearby corners, nearly
parallel rays, and arithmetic exhaustion. These checks are independent
algorithms within the test program; they are not an external agent review
or a sage-flatsurf acceptance receipt.

`parameter` is exact: end = start + parameter × incoming. Directions are
not normalized, and their squared length is preserved by reflection.
The event's length and cumulative distance are floating display estimates.
They never affect predicates or subsequent mathematical trajectories. A
distance-driven renderer may reveal part of a segment without feeding its
rounded screen position back into the mathematical ray. Each ray also has
an explicit maximum number of events.

`billiard_polygon_rectangle` constructs exact geometry from representable
binary width/height values. The existing rectangle controls, launch
position, angle, wall dragging, pause behavior and corner reversal continue
to use the existing rectangle engine. The new mathematical polygon engine
has a separate vertex-stop policy. Surface flow through a regular marked
point and termination at a genuine cone singularity belong to the later
unfolding/flow jobs; a polygon vertex alone does not establish cone order.

No Android integration or general rational-polygon completion is claimed
by this module. ICK remains the only permitted compiler for its C and test
C. The user authorized execution of this campaign without Flexible Pipes
on 19 October 2026; that exception does not change compiler or evidence rules.
