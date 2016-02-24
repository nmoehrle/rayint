workspace "rayint"
    configurations {"debug", "release", "profile"}
    flags {"C++11"}

    configuration "release"
        targetdir "build/release"
        optimize "On"

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

        libdirs {"elibs/mve/libs/**/", os.findlib("png"), os.findlib("jpeg"), os.findlib("tiff")}
        links {"mve", "mve_util", "jpeg", "png", "tiff"}
