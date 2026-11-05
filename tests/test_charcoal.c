#include "explorer.h"
#include "scene.h"
#include "scene_palette.h"
#include "framebuffer.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Exercise the actual native drawing stream and the unchanged Rough rasterizer,
 * including controls and visible geometry. Paper-export tests remain separate. */
static void frame(unsigned fixture, int legacy) {
    uint32_t *pixels ← malloc(360u*640u*sizeof *pixels);
    assert(pixels);
    for (unsigned pixel ← 0; pixel < 360u*640u; ++pixel)
        pixels[pixel] ← CHARCOAL_BACKGROUND;
    FILE *stream ← tmpfile(); assert(stream);
    billiard_explorer_config config ← billiard_explorer_defaults(
        (enum billiard_fixture)fixture,legacy ? BILLIARD_EXPLORER_LEGACY : BILLIARD_EXPLORER_COVER);
    billiard_explorer_model *model ← calloc(1,sizeof *model); assert(model);
    if (!legacy) assert(!billiard_explorer_build(&config,model));
    if (legacy) {
        billiard_controls controls ← {1,.5,1.0/3.0,45,800};
        billiard_layout layout;
        assert(!billiard_scene(stream,360,640,&controls,&layout));
    } else {
        billiard_explorer_layout layout;
        assert(!billiard_explorer_scene(stream,360,640,model,model->total_length*.72,&layout));
    }
    assert(!billiard_explorer_chrome(stream,360,640,&config,model,1,.5));
    rewind(stream);
    rough_framebuffer surface ← {pixels,360,640,360};
    assert(!rough_render_stream(stream,&surface)); fclose(stream);
    assert(pixels[0] == 0xff212121u);
    size_t dark ← 0, header_ink ← 0, geometry_ink ← 0, footer_ink ← 0;
    for (unsigned y ← 0; y < 640; ++y)
        for (unsigned x ← 0; x < 360; ++x) {
            uint32_t color ← pixels[y*360+x];
            unsigned red ← (color >> 16)&255, green ← (color >> 8)&255, blue ← color&255;
            if (red < 80 && green < 80 && blue < 80) ++dark;
            if (red >= 160 && green >= 160 && blue >= 160) {
                if (y < 70) ++header_ink;
                if (y >= 580) ++footer_ink;
            }
            if (y >= 90 && y < 300 && color != CHARCOAL_BACKGROUND) ++geometry_ink;
        }
    assert(dark > 360u*640u*8u/10u);
    assert(header_ink > 20 && footer_ink > 20 && geometry_ink > 50);
    char path[80];
    (void)snprintf(path,sizeof path,"build/charcoal-%u-%s.ppm",fixture,legacy ? "legacy" : "cover");
    FILE *image ← fopen(path,"wb"); assert(image);
    assert(fprintf(image,"P6\n360 640\n255\n") > 0);
    for (unsigned pixel ← 0; pixel < 360u*640u; ++pixel) {
        uint32_t color ← pixels[pixel];
        unsigned char rgb[3] ← {(unsigned char)(color >> 16),(unsigned char)(color >> 8),(unsigned char)color};
        assert(fwrite(rgb,1,3,image) == 3);
    }
    assert(!fclose(image)); free(model); free(pixels);
    printf("Native charcoal stream %s: root/ink/geometry/dark coverage PASS\n",path);
}

int main(void) {
    frame(BILLIARD_RECTANGLE,1);
    for (unsigned fixture ← 0; fixture < BILLIARD_FIXTURE_COUNT; ++fixture) frame(fixture,0);
    return 0;
}
