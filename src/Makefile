# Include makefile setup

include Makefile.setup

# Set TARGET and linker flags for Litmus

LTM_TARGET := litmus.x

LIB_LIST := $(OPENMP_LD)

# Select compiler flags
ifneq "$(MAKECMDGOALS)" "release"
  CXXFLAGS := $(CXX_DEB_FLG)
  CPPFLAGS := $(CPP_DEB_FLG)
else
  CXXFLAGS := $(CXX_REL_FLG)
  CPPFLAGS := $(CPP_REL_FLG)
endif 

# Add includes directories
LTM_INCLUDES := ${shell find * -type d -print}

# Update CXX flags & includes path
CXXFLAGS += $(addprefix -I $(CURDIR)/,$(LTM_INCLUDES))
vpath %.hpp $(LTM_INCLUDES)

# Add list of internal libraries
LTM_LIBS_DIR := ${shell find * -type d -print}

# Generate all internal libraries related files
LTM_LIBS = $(foreach d,$(LTM_LIBS_DIR),$(addsuffix .a,$(addprefix $d/,$d)))
LTM_LIBS_OBJS = $(foreach d,$(LTM_LIBS_DIR),$(addprefix $d/,*.o))
LTM_LIBS_DEPS = $(foreach d,$(LTM_LIBS_DIR),$(addprefix $d/,*.d))
LTM_LIBS_OPTS = $(foreach d,$(LTM_LIBS_DIR),$(addprefix $d/,*.optrpt))
LTM_LIBS_LIST = $(foreach d,$(LTM_LIBS_DIR),$(addprefix -L$d -l,$d))

# Debug
debug: $(LTM_LIBS) litmus.o
	@echo ====== Linking debug version of $(LTM_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(LTM_TARGET) \
	$(LTM_LIBS_OBJS) litmus.o $(LIB_LIST)

# Release
release: $(LTM_LIBS) litmus.o
	@echo ====== Linking release version of $(LTM_TARGET)...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(LTM_TARGET) \
	$(LTM_LIBS) litmus.o $(LIB_LIST)

# Clean up
clean:
	@echo ====== Removing temporary files...
	$(RM) $(LTM_LIBS)
	$(RM) $(LTM_LIBS_OBJS)
	$(RM) $(LTM_LIBS_DEPS)
	$(RM) $(LTM_LIBS_OPTS)

# Set internal libraries generation rule
.PHONY: $(LTM_LIBS) $(LTM_LIBS_DIR)
$(LTM_LIBS): $(LTM_LIBS_DIR)
$(LTM_LIBS_DIR):
	make --directory=$@ $(MAKECMDGOALS)
