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

## spdlog
function(tetwild_download_spdlog)
    tetwild_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.1.0
    )
endfunction()
