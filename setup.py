import setuptools
from setuptools import Extension
from setuptools.command.build_ext import build_ext as _build_ext
import glob
import sys
import os
import platform
import subprocess
import re

try:
    import pybind11
except ImportError:
    print("Error: pybind11 is required to build this package.")
    print("Setup `external/pybind11/include`")
    class pybind11:
        @staticmethod
        def get_include(): return "external/pybind11/include"

extension_name = "pcst_fast"
package_name = "pcst-fast"
version = "1.0.0"

MIN_GCC_VERSION = (13, 0)
MIN_CLANG_VERSION = (15, 0)
MIN_MSVC_VERSION = (19, 34)

class build_ext(_build_ext):
    def run(self):
        self.check_compiler_version()
        _build_ext.run(self)

    def get_compiler_version(self, compiler_executable):
        """Attempts to get the version tuple (major, minor) for a compiler."""
        version_output = None
        try:
            version_output = subprocess.check_output(
                [compiler_executable, '--version'], stderr=subprocess.STDOUT, timeout=5
            ).decode(errors='replace')
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
             if sys.platform == 'win32' and 'cl' in os.path.basename(compiler_executable).lower():
                 try:
                     proc = subprocess.Popen([compiler_executable], stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True, errors='replace')
                     stdout, stderr = proc.communicate(timeout=5)
                     version_output = stderr or stdout
                 except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
                     return None, None
             else:
                 return None, None

        if not version_output:
             return None, None

        m = re.search(r' (\d{1,3})\.(\d{1,3})\.?(\d{1,3})?', version_output)
        if m and ('gcc' in version_output.lower() or 'g++' in version_output.lower()):
            return (int(m.group(1)), int(m.group(2))), version_output

        m = re.search(r'clang version (\d{1,3})\.(\d{1,3})\.?(\d{1,3})?', version_output, re.IGNORECASE)
        if m:
            return (int(m.group(1)), int(m.group(2))), version_output

        m = re.search(r'Version[^\d]*(\d{1,3})\.(\d{1,3})\.?(\d{1,3})?', version_output)
        if m and 'Microsoft' in version_output:
            return (int(m.group(1)), int(m.group(2))), version_output

        m = re.search(r'(\d+)\.(\d+)', version_output)
        if m:
             print(f"Warning: Using generic version parsing for output:\n{version_output}")
             return (int(m.group(1)), int(m.group(2))), version_output

        return None, version_output

    def check_compiler_version(self):
        """Checks if the detected C++ compiler version meets minimum requirements."""
        compiler_executable = None
        if hasattr(self, 'compiler') and self.compiler and hasattr(self.compiler, 'compiler_cxx') and self.compiler.compiler_cxx:
            compiler_executable = self.compiler.compiler_cxx[0]
        else:
            compiler_executable = os.environ.get('CXX', None)
            if not compiler_executable:
                 if sys.platform == 'win32':
                     compiler_executable = 'cl'
                 elif sys.platform == 'darwin':
                     compiler_executable = 'clang++'
                 else:
                     compiler_executable = 'g++'

        print(f"--- Checking C++ compiler version for: {compiler_executable} ---")
        compiler_name = os.path.basename(compiler_executable).lower()
        version_tuple, version_output = self.get_compiler_version(compiler_executable)

        if version_tuple is None:
            print(f"Warning: Could not determine version for compiler '{compiler_executable}'. Output:\n{version_output or 'No output'}")
            print("Proceeding, but build might fail if compiler is too old for C++23.")
            return

        major, minor = version_tuple
        print(f"Found compiler version: {major}.{minor}")

        min_version = None
        compiler_type = "Unknown"

        if 'g++' in compiler_name or (version_output and 'gcc' in version_output.lower() and 'clang' not in version_output.lower()):
            min_version = MIN_GCC_VERSION
            compiler_type = "g++"
        elif 'clang++' in compiler_name or (version_output and 'clang' in version_output.lower()):
            min_version = MIN_CLANG_VERSION
            compiler_type = "clang++"
        elif 'cl.exe' in compiler_name or compiler_name == 'cl':
            min_version = MIN_MSVC_VERSION
            compiler_type = "MSVC"
        else:
            if compiler_executable:
                print(f"Warning: Unrecognized compiler type '{compiler_name}'. Assuming compatibility.")
            else:
                print(f"Warning: Compiler executable could not be determined. Cannot check version.")
            return

        if min_version and (major, minor) < min_version:
             raise RuntimeError(
                 f"{compiler_type} version {min_version[0]}.{min_version[1]} or higher is required for C++23 features. "
                 f"Found {major}.{minor}."
             )
        elif min_version:
            print(f"{compiler_type} version check passed (found {major}.{minor}, requires >={min_version[0]}.{min_version[1]}).")

cpp_sources = glob.glob("src/**/*.cc", recursive=True) + \
              glob.glob("bindings/python/*.cc", recursive=True)

if not cpp_sources:
    raise RuntimeError("Could not find C++ source files. Check src/ and bindings/python directories.")
print(f"Found source files: {cpp_sources}")

include_dirs = [
    "include",
    pybind11.get_include(),
]
print(f"Using include directories: {include_dirs}")

extra_compile_args = [
    "-O3",
    "-Wall",
    "-Wextra",
    "-pedantic",
    "-fPIC",
    "-DNDEBUG",
]
extra_link_args = []

if sys.platform == "win32":
    extra_compile_args = [
        "/O2",
        "/W3",
        "/DNDEBUG",
        "/EHsc",
        "/std:c++latest",
        "/bigobj",
        "/utf-8",
    ]
    extra_link_args = []
elif sys.platform == "darwin":
    extra_compile_args.extend([
        "-std=c++23",
        "-mmacosx-version-min=10.15",
    ])
    extra_link_args.extend([
        "-stdlib=libc++",
        "-mmacosx-version-min=10.15",
    ])
else:
    extra_compile_args.extend([
        "-std=c++23",
    ])
    extra_link_args.extend([
         "-pthread",
    ])

print(f"Using extra_compile_args: {extra_compile_args}")
print(f"Using extra_link_args: {extra_link_args}")

ext_modules = [
    Extension(
        extension_name,
        sources=sorted(cpp_sources),
        include_dirs=include_dirs,
        language='c++',
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "Fast implementation of the Prize-Collecting Steiner Forest / Tree algorithm (PCSF / PCST)."

setuptools.setup(
    name=package_name,
    version=version,
    description="A fast C++ implementation of PCSF/PCST with Python bindings",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: C++",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    python_requires='>=3.9',
    setup_requires=[f'pybind11>=2.10'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
