#ifndef BILLIARD_SCENE_PALETTE_H
#define BILLIARD_SCENE_PALETTE_H

#include "../_charcoal/native/presentation/charcoal.h"

/* Native application chrome uses charcoal. Paper drawing streams retain their
 * authored colors; BILLIARD_CHARCOAL also exercises the native palette on host.
 * This selects drawing ink before rasterization, never filters rendered data. */
#if defined(__ANDROID__) || defined(BILLIARD_CHARCOAL)
#define BILLIARD_RGB(role, red, green, blue) \
    (unsigned)(((role) >> 16) & 255u), \
    (unsigned)(((role) >> 8) & 255u), (unsigned)((role) & 255u)
#define BILLIARD_ALPHA(charcoal, paper) (charcoal)
#define BILLIARD_LABEL(alpha) ((alpha) < .5 ? CHARCOAL_SECONDARY : CHARCOAL_TEXT)
#else
#define BILLIARD_RGB(role, red, green, blue) \
    (unsigned)(red), (unsigned)(green), (unsigned)(blue)
#define BILLIARD_ALPHA(charcoal, paper) (paper)
#define BILLIARD_LABEL(alpha) CHARCOAL_TEXT
#endif

/* Blue identifies the trajectory, as in the original paper scene. */
#define BILLIARD_TRAJECTORY 0xff96c9ebu

#endif
