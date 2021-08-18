
UFBXT_TEST(fuzz_files)
#if UFBXT_IMPL
{
	size_t ok = 0;
	size_t i = g_fuzz_step < SIZE_MAX ? g_fuzz_step : 0;
	for (; i < 10000; i++) {
		char name[512];
		char buf[512];
		snprintf(name, sizeof(name), "fuzz_%04zu", i);
		snprintf(buf, sizeof(buf), "%sfuzz/fuzz_%04zu.fbx", data_root, i);

		size_t size;
		void *data = ufbxt_read_file(buf, &size);
		if (!data) break;

		ufbx_error error;

		ufbx_load_opts load_opts = { 0 };
		load_opts.temp_allocator.memory_limit = 0x4000000; // 64MB
		load_opts.result_allocator.memory_limit = 0x4000000; // 64MB

		ufbx_scene *scene = ufbx_load_memory(data, size, &load_opts, &error);
		if (scene) {
			ufbxt_check_scene(scene);
			ok++;
		}


		ufbx_free_scene(scene);

		ufbx_load_opts stream_opts = load_opts;
		ufbxt_init_allocator(&stream_opts.temp_allocator);
		ufbxt_init_allocator(&stream_opts.result_allocator);
		stream_opts.read_buffer_size = 1;
		stream_opts.temp_allocator.huge_threshold = 1;
		stream_opts.result_allocator.huge_threshold = 1;
		ufbx_scene *streamed_scene = ufbx_load_file(buf, &stream_opts, &error);
		if (streamed_scene) {
			ufbxt_check_scene(streamed_scene);
			ufbxt_assert(scene);
		} else {
			ufbxt_assert(!scene);
		}
		ufbx_free_scene(streamed_scene);

		ufbxt_do_fuzz(scene, streamed_scene, name, data, size);

		free(data);
	}

	ufbxt_logf(".. Loaded fuzz files: %zu (%zu non-errors)", i, ok);
}
#endif
