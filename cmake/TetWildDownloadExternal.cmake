################################################################################
include(DownloadProject)

# Shortcut function
function(tetwild_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${TETWILD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${TETWILD_EXTERNAL}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(tetwild_download_libigl)
    tetwild_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        45cfc79fede992ea3923ded9de3c21d1c4faced1
    )
endfunction()

## geogram
function(tetwild_download_geogram)
    tetwild_download_project(geogram
        GIT_REPOSITORY https://github.com/Yixin-Hu/geogram/
        GIT_TAG        b613750341a6cdd31ae8df80ecfc26ac7ca1a6ad
    )
endfunction()

## fmt
function(tetwild_download_fmt)
    tetwild_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        5.2.0
    )
endfunction()

## spdlog
function(tetwild_download_spdlog)
    tetwild_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.1.0
    )
endfunction()

## CLI11
function(tetwild_download_cli11)
    tetwild_download_project(cli11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.6.1.tar.gz
        URL_MD5 48ef97262adb0b47a2f0a7edbda6e2aa
    )
endfunction()
