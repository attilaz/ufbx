@ctype mat4 hmm_mat4

@vs mesh_vs

in vec3 a_position;
in vec3 a_normal;

out vec3 v_normal;

uniform mesh_vertex_ubo {
    mat4 geometry_to_world;
    mat4 normal_to_world;
    mat4 world_to_clip;
};

void main() {
    vec4 v = vec4(a_position, 1.0);
    v = geometry_to_world * v;
    v = world_to_clip * v;
    gl_Position = v;

    vec4 n = vec4(a_normal, 0.0);
    n = normal_to_world * n;
    v_normal = n.xyz;
}
@end

@fs mesh_fs
in vec3 v_normal;

out vec4 o_color;

void main() {
    vec3 normal = normalize(v_normal);
    float f = dot(normal, normalize(vec3(1.0, 1.0, 1.0))) * 0.5 + 0.5;
    o_color = vec4(f, f, f, 1.0);
}
@end

@program mesh mesh_vs mesh_fs
