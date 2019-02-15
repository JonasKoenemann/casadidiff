// Throw informative error message
#define CASADI_THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

// Throw informative error message
#define CASADI_THROW_ERROR_OBJ(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " for node of type " \
  + this->class_name() + " at " + CASADI_WHERE + ":\n" + std::string(WHAT));


  // Throw informative error message
  #define THROW_ERROR(FNAME, WHAT) \
  throw CasadiException("Error in Function::" FNAME " for '" + this->name() + "' "\
    "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
    + string(WHAT));

  // Throw informative error message from constructor
  #define THROW_ERROR_NOOBJ(FNAME, WHAT, CLASS_NAME) \
  throw CasadiException("Error in Function::" FNAME " for '" + name + "' "\
      "[" CLASS_NAME "] at " + CASADI_WHERE + ":\n"\
      + string(WHAT));
