// Throw informative error message
#define CASADI_THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

// Throw informative error message
#define CASADI_THROW_ERROR_OBJ(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " for node of type " \
  + this->class_name() + " at " + CASADI_WHERE + ":\n" + std::string(WHAT));
