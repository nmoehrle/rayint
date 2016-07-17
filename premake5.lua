workspace "rayint"
    configurations {"debug", "release", "profile"}
    flags {"C++11"}
    location "build"

    configuration "release"
        targetdir "build/release"
        buildoptions {"-fopenmp"}
        optimize "On"

    filter {"system:linux"}
        linkoptions {"-pthread"}

    configuration "debug"
        targetdir "build/debug"
        defines {"Symbols"}

    os.execute("git submodule init")
    os.execute("git submodule update")

    project "raycast"
        kind "ConsoleApp"
        language "C++"
        files {"apps/raycast/*.h", "apps/raycast/*.cpp"}
        includedirs {"libs", "elibs/mve/libs"}

        libdirs {"elibs/mve/libs/**/"}
        links {"gomp", "mve", "mve_util", "jpeg", "png", "tiff"}
