# Vina-GPU 2.1 AMD/RDNA4 patches

Apply to a fresh clone of `DeltaGroupNJUPT/Vina-GPU-2.1`, subdir
`AutoDock-Vina-GPU-2.1/`, then build with `build_vina_gpu_amd.sh`.
Needed because upstream targets NVIDIA + older Boost; AMD's OpenCL compiler and
modern MinGW/Boost are stricter.

## 1. main/main.cpp — Boost convenience.hpp removed in modern Boost
```diff
-#include <boost/filesystem/convenience.hpp> // filesystem::basename
+#include <boost/filesystem.hpp> // filesystem::basename (convenience.hpp removed in modern Boost)
```

## 2. lib/main_procedure_cl.cpp — windows.h for Sleep()
After the `#include <stdio.h>` near the top:
```diff
 #include <stdio.h>
+#ifdef WIN32
+#include <windows.h> // Sleep()
+#endif
```

## 3. lib/main_procedure_cl.cpp — runtime kernel cache
Replace the `#ifdef BUILD_KERNEL_FROM_SOURCE ... #else <load> #endif` block around
the two kernel builds with a runtime check: load `Kernel{1,2}_Opt.bin` if present,
else compile-from-source and save. Opening of the build-from-source block:
```diff
-#ifdef BUILD_KERNEL_FROM_SOURCE
+	{
+	FILE* _c1 = fopen("Kernel1_Opt.bin", "rb");
+	FILE* _c2 = fopen("Kernel2_Opt.bin", "rb");
+	bool _have_bins = (_c1 != NULL && _c2 != NULL);
+	if (_c1) fclose(_c1);
+	if (_c2) fclose(_c2);
+	if (!_have_bins) {
 	const std::string default_work_path = ".";
```
Closing (after `SaveProgramToBinary(programs[1], "Kernel2_Opt.bin");`):
```diff
 	SaveProgramToBinary(programs[1], "Kernel2_Opt.bin");
-#else
-	programs[0] = SetupBuildProgramWithBinary(context, devices, (opencl_binary_path + std::string("/Kernel1_Opt.bin")).c_str());
-	programs[1] = SetupBuildProgramWithBinary(context, devices, (opencl_binary_path + std::string("/Kernel2_Opt.bin")).c_str());
-#endif
+	} else {
+		printf("\n\nLoading cached kernels (gfx1201)"); fflush(stdout);
+		programs[0] = SetupBuildProgramWithBinary(context, devices, "Kernel1_Opt.bin");
+		programs[1] = SetupBuildProgramWithBinary(context, devices, "Kernel2_Opt.bin");
+	}
+	}
```

## 4. OpenCL/src/wrapcl.cpp — drop -Werror (AMD compiler rejects trailing null chars)
Both the LARGE_BOX and SMALL_BOX option strings:
```diff
-    std::string option = " -Werror -cl-single-precision-constant ... ";
+    std::string option = " -w -cl-single-precision-constant ... ";
```

## 5. OpenCL/src/wrapcl.cpp — binary file I/O must be binary mode (THE CL_INVALID_BINARY fix)
```diff
-    FILE* fp = fopen(file_name, "w");          // SaveProgramToBinary
+    FILE* fp = fopen(file_name, "wb");
-    FILE* program_handle = fopen(binary_file_name, "r");   // SetupBuildProgramWithBinary
+    FILE* program_handle = fopen(binary_file_name, "rb");
```
Text mode on Windows translates `\n`<->`\r\n`, corrupting the saved CL binary so
reload fails with `CL_INVALID_BINARY`.

## Build defines
`-DAMD_PLATFORM` (select AMD OpenCL platform), `-DOPENCL_2_0`,
`-DCL_TARGET_OPENCL_VERSION=200`, `-DSMALL_BOX`, `-DBUILD_KERNEL_FROM_SOURCE`.
