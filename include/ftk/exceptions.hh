#ifndef __FTK_EXCEPTIONS_HH
#define __FTK_EXCEPTIONS_HH

#include <ftk/error.hh>
#include <exception>
#include <string>
#include <sstream>

namespace ftk {

/**
 * @brief Base exception class for all FTK errors
 *
 * All FTK exceptions derive from this base class, which provides
 * contextual error information including error code, message, and
 * optional file/line information.
 */
class ftk_error : public std::exception {
public:
  ftk_error(int error_code, const std::string& msg = "")
    : error_code_(error_code), message_(msg)
  {
    build_what_string();
  }

  ftk_error(int error_code, const std::string& msg,
            const std::string& file, int line)
    : error_code_(error_code), message_(msg), file_(file), line_(line)
  {
    build_what_string();
  }

  virtual ~ftk_error() noexcept = default;

  const char* what() const noexcept override {
    return what_.c_str();
  }

  int error_code() const noexcept {
    return error_code_;
  }

  const std::string& message() const noexcept {
    return message_;
  }

  const std::string& file() const noexcept {
    return file_;
  }

  int line() const noexcept {
    return line_;
  }

protected:
  void build_what_string() {
    std::ostringstream oss;
    oss << "[FTK ERROR] " << ftk_err2str(error_code_);
    if (!message_.empty()) {
      oss << ": " << message_;
    }
    if (!file_.empty()) {
      oss << " (at " << file_ << ":" << line_ << ")";
    }
    what_ = oss.str();
  }

  int error_code_;
  std::string message_;
  std::string file_;
  int line_ = -1;
  std::string what_;
};

/**
 * @brief Exception for I/O related errors (file operations, network, etc.)
 */
class io_error : public ftk_error {
public:
  io_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  io_error(int error_code, const std::string& msg,
           const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for mesh-related errors
 */
class mesh_error : public ftk_error {
public:
  mesh_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  mesh_error(int error_code, const std::string& msg,
             const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for ndarray-related errors
 */
class ndarray_error : public ftk_error {
public:
  ndarray_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  ndarray_error(int error_code, const std::string& msg,
                const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for runtime errors (computation failures, unsupported operations)
 */
class runtime_error : public ftk_error {
public:
  runtime_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  runtime_error(int error_code, const std::string& msg,
                const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for logic errors (programming errors, invalid arguments)
 */
class logic_error : public ftk_error {
public:
  logic_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  logic_error(int error_code, const std::string& msg,
              const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for missing optional dependencies
 */
class dependency_error : public ftk_error {
public:
  dependency_error(int error_code, const std::string& msg = "")
    : ftk_error(error_code, msg) {}

  dependency_error(int error_code, const std::string& msg,
                   const std::string& file, int line)
    : ftk_error(error_code, msg, file, line) {}
};

/**
 * @brief Exception for not-yet-implemented features
 */
class not_implemented_error : public ftk_error {
public:
  not_implemented_error(const std::string& msg = "")
    : ftk_error(FTK_ERR_NOT_IMPLEMENTED, msg) {}

  not_implemented_error(const std::string& msg,
                        const std::string& file, int line)
    : ftk_error(FTK_ERR_NOT_IMPLEMENTED, msg, file, line) {}
};

// Convenience macros for throwing with file/line information
#define FTK_THROW(exception_type, error_code, msg) \
  throw exception_type(error_code, msg, __FILE__, __LINE__)

#define FTK_THROW_IO_ERROR(error_code, msg) \
  FTK_THROW(ftk::io_error, error_code, msg)

#define FTK_THROW_MESH_ERROR(error_code, msg) \
  FTK_THROW(ftk::mesh_error, error_code, msg)

#define FTK_THROW_NDARRAY_ERROR(error_code, msg) \
  FTK_THROW(ftk::ndarray_error, error_code, msg)

#define FTK_THROW_RUNTIME_ERROR(error_code, msg) \
  FTK_THROW(ftk::runtime_error, error_code, msg)

#define FTK_THROW_LOGIC_ERROR(error_code, msg) \
  FTK_THROW(ftk::logic_error, error_code, msg)

#define FTK_THROW_DEPENDENCY_ERROR(error_code, msg) \
  FTK_THROW(ftk::dependency_error, error_code, msg)

#define FTK_THROW_NOT_IMPLEMENTED(msg) \
  throw ftk::not_implemented_error(msg, __FILE__, __LINE__)

/**
 * @brief Helper function to throw appropriate exception based on error code
 *
 * This function examines the error code and throws the appropriate
 * exception type (io_error, mesh_error, etc.)
 *
 * This is the exception-based replacement for ftk_fatal() in library code.
 */
inline void ftk_throw(int error_code, const std::string& msg = "",
                      const std::string& file = "", int line = -1) {
  // File I/O errors
  if (error_code >= FTK_ERR_FILE_NOT_FOUND && error_code < FTK_ERR_NOT_BUILT_WITH_ADIOS2) {
    if (file.empty())
      throw io_error(error_code, msg);
    else
      throw io_error(error_code, msg, file, line);
  }

  // Missing dependency errors
  if (error_code >= FTK_ERR_NOT_BUILT_WITH_ADIOS2 && error_code < FTK_ERR_NDARRAY_MULTIDIMENSIONAL_COMPONENTS) {
    if (file.empty())
      throw dependency_error(error_code, msg);
    else
      throw dependency_error(error_code, msg, file, line);
  }

  // Ndarray errors
  if (error_code >= FTK_ERR_NDARRAY_MULTIDIMENSIONAL_COMPONENTS && error_code < FTK_ERR_ACCELERATOR_UNSUPPORTED) {
    if (file.empty())
      throw ndarray_error(error_code, msg);
    else
      throw ndarray_error(error_code, msg, file, line);
  }

  // Mesh errors
  if (error_code >= FTK_ERR_MESH_UNSUPPORTED_FORMAT && error_code <= FTK_ERR_MESH_EMPTY) {
    if (file.empty())
      throw mesh_error(error_code, msg);
    else
      throw mesh_error(error_code, msg, file, line);
  }

  // Not implemented
  if (error_code == FTK_ERR_NOT_IMPLEMENTED) {
    if (file.empty())
      throw not_implemented_error(msg);
    else
      throw not_implemented_error(msg, file, line);
  }

  // Generic runtime error for everything else
  if (file.empty())
    throw runtime_error(error_code, msg);
  else
    throw runtime_error(error_code, msg, file, line);
}

/**
 * @brief Replacement for ftk::fatal() that throws exceptions instead of calling exit()
 *
 * For library code, use this instead of ftk_fatal() from error.hh.
 * CLI tools can continue using ftk_fatal() which calls exit().
 */
namespace detail {
  inline void fatal(int error_code, const std::string& msg = "") {
    ftk_throw(error_code, msg);
  }

  inline void fatal(int error_code, const std::string& msg,
                   const char* file, int line) {
    ftk_throw(error_code, msg, file, line);
  }
}

// For drop-in replacement of ftk::fatal() in existing code
using detail::fatal;

} // namespace ftk

#endif // __FTK_EXCEPTIONS_HH
