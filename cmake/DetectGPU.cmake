find_program(_nvidia_smi "nvidia-smi")
# message ("Detecting NVidia GPUs...")
if (_nvidia_smi)
    set(DETECT_GPU_COUNT_NVIDIA_SMI 0)
    # execute nvidia-smi -L to get a short list of GPUs available
    exec_program(${_nvidia_smi} ARGS -L
        OUTPUT_VARIABLE _nvidia_smi_out
        RETURN_VALUE    _nvidia_smi_ret)
    # process the stdout of nvidia-smi
    if (_nvidia_smi_ret EQUAL 0)
        message(${_nvidia_smi_out})
        # convert string with newlines to list of strings
        string(REGEX REPLACE "\n" ";" _nvidia_smi_out "${_nvidia_smi_out}")
        foreach(_line ${_nvidia_smi_out})
            if (_line MATCHES "^GPU [0-9]+:")
                math(EXPR DETECT_GPU_COUNT_NVIDIA_SMI "${DETECT_GPU_COUNT_NVIDIA_SMI}+1")
                # the UUID is not very useful for the user, remove it
                string(REGEX REPLACE " \\(UUID:.*\\)" "" _gpu_info "${_line}")
                if (NOT _gpu_info STREQUAL "")
                    list(APPEND DETECT_GPU_INFO "${_gpu_info}")
                endif()
            endif()
        endforeach()

        # check_num_gpu_info(${DETECT_GPU_COUNT_NVIDIA_SMI} DETECT_GPU_INFO)
        set(DETECT_GPU_COUNT ${DETECT_GPU_COUNT_NVIDIA_SMI})
    endif()
    # message ("Number of detected NVidia GPUs: " ${DETECT_GPU_COUNT})
endif()
