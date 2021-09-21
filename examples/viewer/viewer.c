
#if defined(_MSC_VER)
#pragma warning(disable: 4116) // (nuklear.h) unnamed type definition in parentheses
#endif

#define NK_IMPLEMENTATION
#define NK_INCLUDE_VERTEX_BUFFER_OUTPUT
#define NK_INCLUDE_FONT_BAKING
#define NK_INCLUDE_DEFAULT_FONT
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_STANDARD_BOOL

#define SOKOL_IMPL
#define SOKOL_GLCORE33

#define arraycount(arr) (sizeof(arr) / sizeof(*(arr)))

#include "external/HandmadeMath.h"
#include "external/nuklear.h"
#include "external/sokol_gfx.h"
#include "external/sokol_app.h"
#include "external/sokol_gl.h"
#include "external/sokol_nuklear.h"
#include "external/sokol_glue.h"
#include "../../ufbx.h"
#include "shaders.h"
#include <stdlib.h>

typedef struct {
    sg_buffer vertex_buffer;
    size_t num_indices;
} viewer_mesh;

typedef struct {
    float position[3];
    float normal[3];
} viewer_vertex;

typedef struct {
    sg_pipeline mesh_pipe;
    sgl_pipeline gl_pipe;
    viewer_mesh *meshes;
} viewer_data;

typedef struct {
    uint32_t selected_element;
} ui_data;

ufbx_scene *scene;
struct nk_context *nk;
viewer_data view;
ui_data ui;

hmm_vec3 vec3_to_hmm(ufbx_vec3 v)
{
    hmm_vec3 dst = { (float)v.x, (float)v.y, (float)v.z };
    return dst;
}

hmm_mat4 matrix_to_hmm(ufbx_matrix mat)
{
    hmm_mat4 dst;
    for (size_t col = 0; col < 4; col++) {
        dst.Elements[col][0] = (float)mat.cols[col].x;
        dst.Elements[col][1] = (float)mat.cols[col].y;
        dst.Elements[col][2] = (float)mat.cols[col].z;
        dst.Elements[col][3] = col == 3 ? 1.0f : 0.0f;
    }
    return dst;
}

void setup_mesh(viewer_mesh *view_mesh, ufbx_mesh *mesh)
{
    sg_destroy_buffer(view_mesh->vertex_buffer);

    size_t num_indices = mesh->num_triangles * 3;

    viewer_vertex *verts = (viewer_vertex*)malloc(sizeof(viewer_vertex) * num_indices);

    size_t dst_ix = 0;
    uint32_t triangles[128];
    for (size_t face_ix = 0; face_ix < mesh->num_faces; face_ix++) {
        size_t num_triangles = ufbx_triangulate_face(triangles, 128, mesh, mesh->faces[face_ix]);
        for (size_t tri_ix = 0; tri_ix < num_triangles * 3; tri_ix++) {
            size_t ix = triangles[tri_ix];
			ufbx_vec3 p = ufbx_get_vertex_vec3(&mesh->skinned_position, ix);
			ufbx_vec3 n = ufbx_get_vertex_vec3(&mesh->skinned_normal, ix);
			verts[dst_ix].position[0] = (float)p.x;
			verts[dst_ix].position[1] = (float)p.y;
			verts[dst_ix].position[2] = (float)p.z;
			verts[dst_ix].normal[0] = (float)n.x;
			verts[dst_ix].normal[1] = (float)n.y;
			verts[dst_ix].normal[2] = (float)n.z;
            dst_ix++;
        }
    }

    view_mesh->num_indices = num_indices;

    sg_buffer_desc desc = { 0 };
    desc.data.ptr = verts;
    desc.data.size = num_indices * sizeof(viewer_vertex);
    view_mesh->vertex_buffer = sg_make_buffer(&desc);
    free(verts);

    ui.selected_element = ~0u;
}

void init(void) {
    sg_setup(&(sg_desc){
        .context = sapp_sgcontext()
    });

    sgl_setup(&(sgl_desc_t){ 0 });

    for (size_t i = 0; i < scene->meshes.count; i++) {
        setup_mesh(&view.meshes[i], scene->meshes.data[i]);
    }

    sg_shader shader = sg_make_shader(mesh_shader_desc(sg_query_backend()));

    {
		sg_pipeline_desc desc = { 0 };
		desc.shader = shader;
		desc.layout.attrs[0].format = SG_VERTEXFORMAT_FLOAT3;
		desc.layout.attrs[1].format = SG_VERTEXFORMAT_FLOAT3;
		desc.depth.write_enabled = true;
		desc.depth.compare = SG_COMPAREFUNC_LESS_EQUAL;
		view.mesh_pipe = sg_make_pipeline(&desc);
	}

    {
		sg_pipeline_desc desc = { 0 };
		desc.depth.write_enabled = true;
		desc.depth.compare = SG_COMPAREFUNC_LESS_EQUAL;
		view.gl_pipe = sgl_make_pipeline(&desc);
	}

    snk_setup(&(snk_desc_t){ 0 });
}

void ui_node_tree(ufbx_node *node)
{
    bool was_selected = ui.selected_element == node->element.id;
    bool is_selected = was_selected;
    const char *name = node->name.data;
    if (!name[0]) name = node->parent ? "(empty name)" : "(root)";
	if (nk_tree_element_push_id(nk, NK_TREE_NODE, name, nk_false, &is_selected, (int)node->element.id)) {
        for (size_t child_ix = 0; child_ix < node->children.count; child_ix++) {
            ufbx_node *child = node->children.data[child_ix];
            ui_node_tree(child);
        }
        nk_tree_element_pop(nk);
    }
	if (is_selected && !was_selected) {
        ui.selected_element = node->element.id;
	}
}

void gl_draw_joint(hmm_vec3 src, hmm_vec3 dst, hmm_vec3 up, float radius, bool selected)
{
    hmm_vec3 axis_z = HMM_SubtractVec3(dst, src);
    float dist = HMM_LengthVec3(axis_z);
    if (dist < 0.00001f) return;
    axis_z = HMM_DivideVec3f(axis_z, dist);
    hmm_vec3 axis_x = HMM_NormalizeVec3(HMM_Cross(axis_z, up));
    hmm_vec3 axis_y = HMM_Cross(axis_z, axis_x);
    radius = HMM_MIN(radius, dist * 0.3f);

    hmm_vec2 silhouette[] = {
        { 0.0f, 0.0f },
        { radius, radius },
        { 0.0f, dist },
    };

    hmm_vec3 axes[] = {
        axis_x,
        axis_y,
        HMM_MultiplyVec3f(axis_x, -1.0f),
        HMM_MultiplyVec3f(axis_y, -1.0f),
    };

    sgl_load_pipeline(view.gl_pipe);

	sgl_begin_quads();

    for (size_t axis_ix = 0; axis_ix < arraycount(axes); axis_ix++) {
        hmm_vec3 axis0 = axes[axis_ix];
        hmm_vec3 axis1 = axes[(axis_ix + 1) % arraycount(axes)];
		for (size_t sil_ix = 0; sil_ix < arraycount(silhouette) - 1; sil_ix++) {
			hmm_vec2 sil0 = silhouette[sil_ix];
			hmm_vec2 sil1 = silhouette[sil_ix + 1];

            hmm_vec3 z0 = HMM_AddVec3(src, HMM_MultiplyVec3f(axis_z, sil0.Y));
            hmm_vec3 z1 = HMM_AddVec3(src, HMM_MultiplyVec3f(axis_z, sil1.Y));
            hmm_vec3 p00 = HMM_AddVec3(z0, HMM_MultiplyVec3f(axis0, sil0.X));
            hmm_vec3 p01 = HMM_AddVec3(z0, HMM_MultiplyVec3f(axis1, sil0.X));
            hmm_vec3 p10 = HMM_AddVec3(z1, HMM_MultiplyVec3f(axis0, sil1.X));
            hmm_vec3 p11 = HMM_AddVec3(z1, HMM_MultiplyVec3f(axis1, sil1.X));

            hmm_vec3 n = HMM_Cross(HMM_SubtractVec3(p11, p00), HMM_SubtractVec3(p10, p01));
            float v = HMM_DotVec3(HMM_NormalizeVec3(n), HMM_Vec3(0.57735f, 0.57735f, 0.57735f)) * 0.3f + 0.5f;
            hmm_vec3 color = HMM_MultiplyVec3f(HMM_Vec3(0.8f, 0.6f, 0.4f), v);
            if (selected) {
                color = HMM_Vec3(0.9f, 0.3f, 0.3f);
            }

            sgl_c3f(color.X, color.Y, color.Z);
            sgl_v3f(p00.X, p00.Y, p00.Z);
            sgl_v3f(p01.X, p01.Y, p01.Z);
            sgl_v3f(p11.X, p11.Y, p11.Z);
            sgl_v3f(p10.X, p10.Y, p10.Z);
        }
    }

    sgl_end();
}

void frame(void) {
    nk = snk_new_frame();

    if (nk_begin(nk, "Test", nk_rect(1000.0f, 0.0f, 400.0f, 1000.0f), 0)) {
        ui_node_tree(scene->root_node);
    }
    nk_end(nk);

    {
		sg_pass_action pass_action = {
			.colors[0] = { .action=SG_ACTION_CLEAR, .value={0.1f, 0.1f, 0.1f, 1.0f} }
		};
		sg_begin_default_pass(&pass_action, sapp_width(), sapp_height());
	}

    hmm_mat4 view_to_clip = HMM_Perspective(80.0f, 1.0f, 0.01f, 100.0f);
    hmm_mat4 world_to_view = HMM_LookAt(HMM_Vec3(1.0f, 1.5f, 1.0f), HMM_Vec3(0.0f, 0.5f, 0.0f), HMM_Vec3(0.0f, 1.0f, 0.0f));
    hmm_mat4 world_to_clip = HMM_MultiplyMat4(view_to_clip, world_to_view);

    sgl_matrix_mode_modelview();
    sgl_load_matrix(&world_to_view.Elements[0][0]);
    sgl_matrix_mode_projection();
    sgl_load_matrix(&view_to_clip.Elements[0][0]);

    for (size_t mesh_ix = 0; mesh_ix < scene->meshes.count; mesh_ix++) {
        ufbx_mesh *mesh = scene->meshes.data[mesh_ix];
        viewer_mesh *view_mesh = &view.meshes[mesh_ix];

        for (size_t inst_ix = 0; inst_ix < mesh->instances.count; inst_ix++) {
            ufbx_node *node = mesh->instances.data[inst_ix];

            sg_apply_pipeline(view.mesh_pipe);

            mesh_vertex_ubo_t vu;
            if (mesh->skinned_is_local) {
                vu.geometry_to_world = matrix_to_hmm(node->geometry_to_world);
                vu.normal_to_world = matrix_to_hmm(ufbx_matrix_for_normals(&node->geometry_to_world));
            } else {
                vu.geometry_to_world = matrix_to_hmm(ufbx_identity_matrix);
                vu.normal_to_world = matrix_to_hmm(ufbx_identity_matrix);
            }
            vu.world_to_clip = world_to_clip;
            sg_apply_uniforms(SG_SHADERSTAGE_VS, SLOT_mesh_vertex_ubo, SG_RANGE_REF(vu));

            sg_bindings binds = {
                .vertex_buffers[0] = view_mesh->vertex_buffer,
            };
            sg_apply_bindings(&binds);

            sg_draw(0, (int)view_mesh->num_indices, 1);
        }
    }

#if 0
    for (size_t bone_ix = 0; bone_ix < scene->bones.count; bone_ix++) {
        ufbx_bone *bone = scene->bones.data[bone_ix];
        for (size_t inst_ix = 0; inst_ix < bone->instances.count; inst_ix++) {
            ufbx_node *node = bone->instances.data[inst_ix];
            ufbx_node *parent = node->parent;
            if (!parent || !parent->bone) continue;

            hmm_vec3 src = vec3_to_hmm(parent->node_to_world.cols[3]);
            hmm_vec3 dst = vec3_to_hmm(node->node_to_world.cols[3]);
            hmm_vec3 up = HMM_NormalizeVec3(vec3_to_hmm(parent->node_to_world.cols[1]));

            bool selected = node->element.id == ui.selected_element;
            gl_draw_joint(src, dst, up, (float)parent->bone->radius, selected);
        }
    }
#endif

    {
        sg_end_pass();
		sg_pass_action pass_action = {
            .colors[0] = { .action=SG_ACTION_LOAD },
            .depth = { .action=SG_ACTION_CLEAR, .value=1.0f },
		};
		sg_begin_default_pass(&pass_action, sapp_width(), sapp_height());
	}

    sgl_draw();

    snk_render(sapp_width(), sapp_height());

    sg_end_pass();
    sg_commit();
}

void event(const sapp_event *ev) {
    snk_handle_event(ev);
}

void cleanup(void) {
    sg_shutdown();
}

sapp_desc sokol_main(int argc, char* argv[]) {
    if (argc < 2) {
        exit(1);
    }

    ufbx_load_opts opts = { 0 };
    opts.evaluate_skinning = true;

    scene = ufbx_load_file(argv[1], &opts, NULL);
    if (!scene) {
        exit(1);
    }
    view.meshes = (viewer_mesh*)calloc(sizeof(viewer_mesh), scene->meshes.count);

    return (sapp_desc){
        .init_cb = init,
        .frame_cb = frame,
		.event_cb = event,
        .cleanup_cb = cleanup,
        .width = 1400,
        .height = 1000,
        .window_title = "Clear Sample",
    };
}
