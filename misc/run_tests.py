#!/usr/bin/env python3

import asyncio
import itertools
import subprocess
import time
import re
import os
import sys
import shutil
import functools
import platform

color_out = sys.stdout

has_color = False
if sys.platform == "win32":
    try:
        import colorama
        colorama.init(wrap=False)
        color_out = colorama.AnsiToWin32(color_out).stream
        has_color = True
    except ImportError:
        print("# Run 'pip install colorama' for colorful output")
else:
    has_color = sys.stdout.isatty()

STYLE_FAIL = "\x1b[31m"
STYLE_WARN = "\x1b[33m"
STYLE_CMD = "\x1b[36m"

def log(line, *, style=""):
    if style and has_color:
        line = style + line + "\x1b[39m"
    print(line, file=color_out, flush=True)

def log_cmd(line):
    log(line, style=STYLE_CMD)

def log_mkdir(path):
    log_cmd("mkdir " + path)

def log_comment(line, fail=False, warn=False):
    style = ""
    if fail: style = STYLE_FAIL
    if warn: style = STYLE_WARN
    if sys.platform == "win32":
        log("rem " + line, style=style)
    else:
        log("# " + line, style=style)

if sys.platform == "win32":
    loop = asyncio.ProactorEventLoop()
else:
    loop = asyncio.get_event_loop()

def flatten_str_list(str_list):
    """Flatten arbitrarily nested str list `item` to `dst`"""
    def inner(result, str_list):
        if isinstance(str_list, str):
            result.append(str_list)
        else:
            for s in str_list:
                inner(result, s)
    result = []
    inner(result, str_list)
    return result

if sys.version_info < (3,8):
    cmd_sema = asyncio.Semaphore(os.cpu_count(), loop=loop)
else:
    cmd_sema = asyncio.Semaphore(os.cpu_count())

async def run_cmd(*args, realtime_output=False, env=None):
    """Asynchronously run a command"""

    await cmd_sema.acquire()

    cmd_args = flatten_str_list(args)
    cmd = cmd_args[0]
    cmd_args = cmd_args[1:]

    pipe = None if realtime_output else asyncio.subprocess.PIPE
    cmdline = subprocess.list2cmdline([cmd] + cmd_args)

    out = err = ""
    ok = False

    log_cmd(cmdline)

    begin = time.time()

    try:
        proc = await asyncio.create_subprocess_exec(cmd, *cmd_args,
            stdout=pipe, stderr=pipe, env=env)

        if not realtime_output:
            out, err = await proc.communicate()
            out = out.decode("utf-8", errors="ignore").strip()
            err = err.decode("utf-8", errors="ignore").strip()

        ok = proc.returncode == 0
    except FileNotFoundError:
        err = f"{cmd} not found"
    except OSError as e:
        err = str(e)

    end = time.time()

    cmd_sema.release()

    return ok, out, err, cmdline, end - begin

class Compiler:
    def __init__(self, name, exe):
        self.name = name
        self.exe = exe
        self.env = None
        self.compile_archs = set()
        self.run_archs = set()

    def run(self, *args, **kwargs):
        return run_cmd(self.exe, args, env=self.env, **kwargs)

class CLCompiler(Compiler):
    def __init__(self, name, exe):
        super().__init__(name, exe)
        self.has_c = True
        self.has_cpp = True

    async def check_version(self):
        _, out, err, _, _ = await self.run()
        m = re.search(r"Version ([.0-9]+) for (\w+)", out + err, re.M)
        if not m: return False
        self.arch = m.group(2).lower()
        self.version = m.group(1)
        return True
    
    def supported_archs(self):
        if self.arch == "x86":
            return ["x86"]
        if self.arch == "x64":
            return ["x64"]
        return []
    
    def compile(self, config):
        sources = config["sources"]
        output = config["output"]

        args = []
        args += sources
        args += ["/MT", "/nologo"]

        obj_dir = os.path.dirname(output)
        args.append(f"/Fo{obj_dir}\\")

        if config.get("warnings", False):
            args.append("/W4")
            args.append("/WX")
        else:
            args.append("/W3")

        if config.get("optimize", False):
            args.append("/Ox")

        if config.get("openmp", False):
            args.append("/openmp")

        args.append("/link")

        args += ["/opt:ref"]
        args.append(f"-out:{output}")

        return self.run(args)

class GCCCompiler(Compiler):
    def __init__(self, name, exe, cpp, has_m32=True):
        super().__init__(name, exe)
        self.has_c = not cpp
        self.has_cpp = cpp
        self.has_m32 = has_m32

    async def check_version(self):
        _, vout, _, _, _ = await self.run("-dumpversion")
        _, mout, _, _, _ = await self.run("-dumpmachine")
        if not (vout and mout): return False
        self.version = vout
        self.arch = mout.lower()
        return True

    def supported_archs(self):
        if "x86_64" in self.arch:
            return ["x86", "x64"] if self.has_m32 else ["x64"]
        if "i686" in self.arch:
            return ["x86"]
        return []
    
    def compile(self, config):
        sources = config["sources"]
        output = config["output"]

        args = []

        if config.get("warnings", False):
            args += ["-Wall", "-Wextra", "-Werror"]

        if config.get("optimize", False):
            args.append("-O2")

        if config.get("openmp", False):
            args.append("-openmp")
        
        if self.has_m32 and config.get("arch", "") == "x86":
            args.append("-m32")

        if self.has_cpp:
            args.append("-std=c++11")
        else:
            args.append("-std=gnu99")
        
        args += sources

        if "msvc" not in self.arch:
            args.append("-lm")

        args += ["-o", output]

        return self.run(args)

class ClangCompiler(GCCCompiler):
    def __init__(self, name, exe, cpp, **kwargs):
        super().__init__(name, exe, cpp, **kwargs)

class EmscriptenCompiler(ClangCompiler):
    def __init__(self, name, exe, cpp):
        super().__init__(name, exe, cpp)

class ZigCompiler(ClangCompiler):
    def __init__(self, name, exe, cpp):
        super().__init__(name, exe, cpp)

@functools.lru_cache(8)
def get_vcvars(bat_name):
    vswhere_path = r"%ProgramFiles(x86)%/Microsoft Visual Studio/Installer/vswhere.exe"
    vswhere_path = os.path.expandvars(vswhere_path)
    if not os.path.exists(vswhere_path):
        raise EnvironmentError("vswhere.exe not found at: %s", vswhere_path)

    vs_path = os.popen('"{}" -latest -property installationPath'.format(vswhere_path)).read().rstrip()
    vsvars_path = os.path.join(vs_path, f"VC\\Auxiliary\\Build\\{bat_name}")

    output = os.popen('"{}" && set'.format(vsvars_path)).read()
    env = { }
    for line in output.splitlines():
        items = tuple(line.split("=", 1))
        if len(items) == 2:
            env[items[0]] = items[1]
    return env

class VsCompiler(Compiler):
    def __init__(self, name, bat, inner):
        super().__init__(name, inner.exe)
        self.bat = bat
        self.inner = inner

    async def check_version(self):
        try:
            env = get_vcvars(self.bat)
            self.inner.env = env
            for key, value in env.items():
                if key.upper() != "PATH": continue
                for path in value.split(";"):
                    cl_path = os.path.join(path, self.exe)
                    if os.path.exists(cl_path):
                        self.inner.exe = cl_path
                        if await self.inner.check_version():
                            self.exe = self.inner.exe
                            self.arch = self.inner.arch
                            self.version = self.inner.version
                            self.has_c = self.inner.has_c
                            self.has_cpp = self.inner.has_cpp
                            return True
            return False
        except:
            return False
    
    def supported_archs(self):
        return self.inner.supported_archs()
    
    def compile(self, config):
        return self.inner.compile(config)

all_compilers = [
    CLCompiler("cl", "cl.exe"),
    GCCCompiler("gcc", "gcc", False),
    GCCCompiler("gcc", "g++", True),
    ClangCompiler("clang", "clang", False),
    ClangCompiler("clang", "clang++", True),
    VsCompiler("vs_cl64", "vcvars64.bat", CLCompiler("cl", "cl.exe")),
    VsCompiler("vs_cl32", "vcvars32.bat", CLCompiler("cl", "cl.exe")),
    # VsCompiler("vs_clang64", "vcvars64.bat", ClangCompiler("clang", "clang.exe", False, has_m32=False)),
    # VsCompiler("vs_clang64", "vcvars64.bat", ClangCompiler("clang", "clang++.exe", True, has_m32=False)),
    # VsCompiler("vs_clang32", "vcvars32.bat", ClangCompiler("clang", "clang.exe", False, has_m32=False)),
    # VsCompiler("vs_clang32", "vcvars32.bat", ClangCompiler("clang", "clang++.exe", True, has_m32=False)),
    # EmscriptenCompiler("emcc", emcc", False),
    # EmscriptenCompiler("emcc", emcc++", True),
]

def gather(aws):
    return asyncio.gather(*aws)

ichain = itertools.chain.from_iterable

async def check_compiler(compiler):
    if await compiler.check_version():
        return [compiler]
    else:
        return []

async def find_compilers():
    return list(ichain(await gather(check_compiler(c) for c in all_compilers)))

class Target:
    def __init__(self, name, compiler, config):
        self.name = name
        self.compiler = compiler
        self.config = config
        self.skipped = False
        self.compiled = False
        self.ran = False
        self.ok = True
        self.log = []

async def compile_target(t):
    arch_test = t.config.get("arch_test", False)
    if t.config["arch"] not in t.compiler.compile_archs and not arch_test:
        t.skipped = True
        return

    path = os.path.dirname(t.config["output"])
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        log_mkdir(path)

    ok, out, err, cmdline, time = await t.compiler.compile(t.config)

    t.log.append("$ " + cmdline)
    t.log.append(out)
    t.log.append(err)

    if ok:
        t.compiled = True
    else:
        t.ok = False

    head = f"Compile {t.name}"
    tail = f"[{time:.1f}s OK]" if ok else ("[WARN]" if arch_test else "[FAIL]")
    log_comment(f"{tail} {head}",
        fail=not ok and not arch_test,
        warn=not ok and arch_test)
    return t

async def run_target(t, args):
    if not t.compiled: return
    arch_test = t.config.get("arch_test", False)
    if t.config["arch"] not in t.compiler.run_archs and not arch_test:
        return

    ok, out, err, cmdline, time = await run_cmd(t.config["output"], args)

    t.log.append("$ " + cmdline)
    t.log.append(out)
    t.log.append(err)

    if ok:
        t.ran = True
    else:
        t.ok = False

    head = f"Run {t.name}"
    tail = f"[{time:.1f}s OK]" if ok else ("[WARN]" if arch_test else "[FAIL]")
    log_comment(f"{tail} {head}",
        fail=not ok and not arch_test,
        warn=not ok and arch_test)
    return t

async def compile_and_run_target(t, args):
    await compile_target(t)
    await run_target(t, args)
    return t

def copy_file(src, dst):
    shutil.copy(src, dst)
    if sys.platform == "win32":
        log_cmd(f"copy {src} {dst}")
    else:
        log_cmd(f"cp {src} {dst}")

exit_code = 0

def decorate_arch(compiler, arch):
    if arch not in compiler.compile_archs:
        return arch + " (FAIL)"
    if arch not in compiler.run_archs:
        return arch + " (compile only)"
    return arch

async def main():
    global exit_code

    log_comment("-- Searching for compilers --")
    compilers = await find_compilers()

    all_configs = {
        "optimize": {
            "debug": { "optimize": False },
            "release": { "optimize": True },
        },
        "arch": {
            "x86": { "arch": "x86" },
            "x64": { "arch": "x64" },
        },
    }

    arch_configs = { "arch": all_configs["arch"] }

    build_path = "build"
    if not os.path.exists(build_path):
        os.makedirs(build_path, exist_ok=True)
        log_mkdir(build_path)

    def compile_permutations(prefix, config, config_options, run_args):
        opt_combos = [[(name, opt) for opt in opts] for name, opts in config_options.items()]
        opt_combos = list(itertools.product(*opt_combos))

        is_cpp = config.get("cpp", False)
        is_c = not is_cpp

        for compiler in compilers:
            if is_c and not compiler.has_c: continue
            if is_cpp and not compiler.has_cpp: continue

            for opts in opt_combos:
                optstring = "_".join(opt for _,opt in opts)
                name = f"{prefix}_{compiler.name}_{optstring}"

                path = os.path.join(build_path, name)

                conf = dict(config)
                conf["output"] = os.path.join(path, config.get("output", "a.exe"))
                for opt_name, opt in opts:
                    conf.update(config_options[opt_name][opt])

                target = Target(name, compiler, conf)
                yield compile_and_run_target(target, run_args)

    exe_suffix = ""
    if sys.platform == "win32":
        exe_suffix = ".exe"

    ctest_tasks = []

    ctest_config = {
        "sources": ["misc/compiler_test.c"],
        "output": "ctest" + exe_suffix,
        "arch_test": True,
    }
    ctest_tasks += compile_permutations("ctest", ctest_config, arch_configs, ["1.5"])

    cpptest_config = {
        "sources": ["misc/compiler_test.cpp"],
        "output": "cpptest" + exe_suffix,
        "cpp": True,
        "arch_test": True,
    }
    ctest_tasks += compile_permutations("cpptest", cpptest_config, arch_configs, ["1.5"])

    bad_compiles = 0
    bad_runs = 0

    compiler_test_tasks = await gather(ctest_tasks)
    for target in compiler_test_tasks:
        arch = target.config["arch"]
        if target.compiled:
            target.compiler.compile_archs.add(arch)
        else:
            bad_compiles += 1
        if target.ran and "sin(1.50) = 1.00" in target.log:
            target.compiler.run_archs.add(arch)
        else:
            bad_runs += 1

    log_comment("-- Detected the following compilers --")

    for compiler in compilers:
        archs = ", ".join(decorate_arch(compiler, arch) for arch in compiler.supported_archs())
        log_comment(f"  {compiler.exe}: {compiler.arch} {compiler.version} [{archs}]")

    log_comment("-- Compiling and running tests --")

    target_tasks = []

    runner_config = {
        "sources": ["test/runner.c", "ufbx.c"],
        "output": "runner" + exe_suffix,
    }
    target_tasks += compile_permutations("runner", runner_config, all_configs, ["-d", "data"])

    cpp_config = {
        "sources": ["misc/test_build.cpp"],
        "output": "cpp" + exe_suffix,
        "cpp": True,
        "warnings": True,
    }
    target_tasks += compile_permutations("cpp", cpp_config, arch_configs, [])

    targets = await gather(target_tasks)

    num_fail = 0
    for target in targets:
        if target.ok: continue
        print()
        log(f"-- FAIL: {target.name}", style=STYLE_FAIL)
        print("\n".join(target.log))
        num_fail += 1

    print()
    if bad_compiles > 0:
        print(f"WARNING: {bad_compiles} pre-test configurations failed to compile")
    if bad_runs > 0:
        print(f"WARNING: {bad_runs} pre-test configurations failed to run")
    print(f"{len(targets) - num_fail}/{len(targets)} targets succeeded")
    if num_fail > 0:
        exit_code = 1

if __name__ == "__main__":
    loop.run_until_complete(main())
    loop.close()
    sys.exit(exit_code)