# Finite translation covers — BIL-UNFOLD-02

Reflecting the polygon instead of its moving point lets the point continue
straight. A copy records an orthogonal matrix from the finite group generated
by reflections parallel to the polygon sides. Across a side, the next matrix
is the old matrix followed by that side reflection. Corresponding endpoints
then differ by one translation. That translation is the chart transition;
the flow vector will not reflect in the cover.

The code enumerates this group exactly in the polygon’s quadratic field, with at most 24 copies.
The number of copies must equal twice the least common multiple of the reduced
angle denominators. Unsupported or uncertified polygons cannot manufacture a
finite cover. Arithmetic and copy limits leave the caller's output unchanged.

All side labels refer to the original polygon. A negative-determinant copy
traverses that boundary in reverse when viewed with the surface orientation.
Glued copies have opposite determinant signs, so their oriented boundaries
are paired in reverse even though corresponding original vertex labels match.
The validator checks both endpoints, inverse pairings, matrix relations,
orthogonality, distinct copies, orientation, and connectedness.

The original billiard, its orientation double and this translation cover are
different objects. The square cover here has four unit-square charts, area 4,
and torus periods 2 and 2i. It is not the one-chart unit-area square torus.
Likewise, the L billiard yields four L charts with area 12, genus two and one
6π cone point. It is not silently replaced by the area-3 three-square surface.
Regular vertex classes are retained as marked points of cone angle 2π and
zero order zero. They are not singularities.

Corner identifications give V. Each paired side gives an edge, and each copy
gives a face. Euler characteristic V − E + F determines genus. Summing the
exact polygon angles at each identified corner gives its cone angle and zero
order; the independent identity Σ zero orders = 2g − 2 checks the result.
The test program reconstructs corner cycles by graph traversal independently
of the production union/find calculation. Hand-derived controls are recorded
in `tests/surface-controls.tsv`; malformed pairings and metadata must fail.

The TSV serializer records the model, all seam translations, cone classes,
genus and exact area. A drawing in one plane is not a global chart or a
topology certificate. Sage/flatsurf is not an application dependency, and no
external Sage run or separate-worker acceptance is claimed by these tests.

## Eighth-angle polygons

The polygon field also selects the chart arithmetic. For the triangle
(0,0), (1,0), (1,√2−1), the three generating reflections act on a
dihedral chart (k,s) by (k+s·j mod 8, −s), where j = 0, 4, 1 for the
horizontal, vertical and sloping sides. The 16 charts yield 24 glued
edges and six vertex classes. Four classes lie over the right-angle
corner. The 3π/8 corner has one lift with total angle 6π, contributing
the single order-two zero of a genus-two surface. The exact total area
is 8√2−8. The test independently checks all 48 directed seam transitions.

Triangle corners can contribute odd numbers of π/8. Topology sums in
units of π/24 before emitting the existing whole-cone angle in units of
π/12. The Q(√2) surface report starts with `radicand<TAB>2`. The
existing period and cylinder routines continue to require Q(√3)
charts.
