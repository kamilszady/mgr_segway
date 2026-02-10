# Copyright 1994-2015 The MathWorks, Inc. / INTECO LTD 2020
#
# File    : UnTransZynq.tmf   v.1.1 2020-10-25 (2K)
#
# Abstract:
#	Template makefile used to create <model>.mak file needed for Linaro.
#
# 	Note that this template is automatically customized by the build 
#       procedure to create "<model>.mk" which then is passed to nmake to 
#       produce <model>.mak
#
#       The following defines can be used to modify the behavior of the
#	build:
#    OPTS           - Additional user defines.
#    USER_SRCS      - Additional user sources, such as files needed by
#                     S-functions.
#    USER_INCLUDES  - Additional include paths 
#                     (i.e. USER_INCLUDES="/I where-ever1 /I where-ever2")

#------------------------ Macros read by make_rtw -----------------------------
#
# The following macros are read by the build procedure:
#
#  MAKECMD         - This is the command used to invoke the make utility
#  HOST            - What platform this template makefile is targeted for 
#                    (i.e. PC or UNIX)
#  BUILD           - Invoke make from the build procedure (yes/no)?
#  SYS_TARGET_FILE - Name of system target file.
#  RTDACZYNQ_TOOLBOX_PATH - path to the toolbox files
#  ELF_NAME        - executable file name


MAKECMD                = D:\Matlab\toolbox\UnTransZynq\tools\make
HOST                   = PC
BUILD                  = yes
SYS_TARGET_FILE        = UnTransZynq.tlc
RTDACZYNQ_TOOLBOX_PATH = D:\Matlab\toolbox\UnTransZynq
ELF_NAME               = uZSimulink.elf


#---------------------- Tokens expanded by make_rtw ---------------------------
#
# The following tokens, when wrapped with "|>" and "<|" are expanded by the 
# build procedure.
#   
#  MODEL_NAME          - Name of the Simulink block diagram
#  MODEL_MODULES       - Any additional generated source modules
#  MAKEFILE_NAME       - Name of makefile created from template makefile <model>.mk
#  MATLAB_ROOT         - Path to where MATLAB is installed. 
#  MATLAB_BIN          - Path to MATLAB executable.
#  SOLVER              - Solver source file name
#  NUMST               - Number of sample times
#  TID01EQ             - yes (1) or no (0): Are sampling rates of continuous task 
#                        (tid=0) and 1st discrete task equal.
#  NCSTATES            - Number of continuous states
#  BUILDARGS           - Options passed in at the command line.
#  MULTITASKING        - yes (1) or no (0): Is solver mode multitasking
#  EXT_MODE            - yes (1) or no (0): Build for external mode
#  EXTMODE_TRANSPORT   - Index of transport mechanism (e.g. tcpip, serial) for extmode
#  EXTMODE_STATIC      - yes (1) or no (0): Use static instead of dynamic mem alloc.
#  EXTMODE_STATIC_SIZE - Size of static memory allocation buffer.

#######################################################
# 
MODEL               = UnTrans_basicLQR
MODULES             = rt_matrx.c rt_printf.c rt_logging.c UnTrans_basicLQR_data.c rtGetInf.c rtGetNaN.c rt_nonfinite.c rtiostream_utils.c microZedLED.c microZedStats.c microZedSwitch.c microZedTimer.c microZedUnTrans.c
MAKEFILE            = UnTrans_basicLQR.mk
MATLAB_ROOT         = D:\Matlab
ALT_MATLAB_ROOT     = D:\Matlab
MATLAB_BIN          = D:\Matlab\bin
ALT_MATLAB_BIN      = D:\Matlab\bin
MASTER_ANCHOR_DIR   = 
START_DIR           = C:\Users\Kszad\Desktop\KAMILS~1\STUDIA~1\air2sem2\mgr\PROJEK~1\testy0
SOLVER              = 
NUMST               = 1
TID01EQ             = 0
NCSTATES            = 0
BUILDARGS           =  TMW_EXTMODE_TESTING=0 MAT_FILE=1 COMBINE_OUTPUT_UPDATE_FCNS=1 INCLUDE_MDL_TERMINATE_FCN=1 MULTI_INSTANCE_CODE=0 MODELREF_TARGET_TYPE=NONE RELATIVE_PATH_TO_ANCHOR=.. OPTS="-DEXT_MODE -DON_TARGET_WAIT_FOR_START=0 -DTID01EQ=0"
CLASSIC_INTERFACE   = 0 
MULTITASKING        = 0
EXT_MODE            = 1
EXTMODE_TRANSPORT   = 0
EXTMODE_STATIC      = 0
EXTMODE_STATIC_SIZE = 1000000
TARGET_LANG_EXT     = c

WAIT_FOR_START_PACKAGE = 1

#---------------------------- Tool Definitions -------------------------------
CC      = $(RTDACZYNQ_TOOLBOX_PATH)\tools\linaroeclipse\bin\arm-xilinx-eabi-gcc 
LD      = $(RTDACZYNQ_TOOLBOX_PATH)\tools\linaroeclipse\bin\arm-xilinx-eabi-gcc 
AR      = $(RTDACZYNQ_TOOLBOX_PATH)\tools\linaroeclipse\bin\arm-xilinx-eabi-ar 

# makefile name - required to create dependencied on -D compile options
MAKEFILE_NAME = $(MODEL)_Times.mk

#------------------------------ Include Path ----------------------------------
COMPILER_INCLUDES = 
INCLUDES = $(MATLAB_INCLUDES) $(COMPILER_INCLUDES) $(USER_INCLUDES) 

#-------------------------- Source Path ---------------------------------------
# User source path
ifdef USER_PATH
EXTRA_PATH = ";$(USER_PATH)"
else
EXTRA_PATH = 
endif

# Path to Toolbox S-Functions
TOOLBOXES  = $(MATLAB_ROOT)\toolbox\comm\commsfun;
TOOLBOXES += $(MATLAB_ROOT)\toolbox\dspblks\dspmex;
TOOLBOXES += $(MATLAB_ROOT)\toolbox\fixpoint;
TOOLBOXES += $(MATLAB_ROOT)\toolbox\fuzzy\fuzzy\src

#--------------------------------- Rules --------------------------------------
include $(MODEL)_Times.mk


LDFLAGS := -Wl,-T -Wl,$(RTDACZYNQ_TOOLBOX_PATH)\src\lscript.ld 
LDFLAGS += -L$(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\lib
LDFLAGS += -L.

CFLAGS  = -I..\$(MODEL)_uZ_rtw
CFLAGS += -I$(MATLAB_ROOT)\extern\include  
CFLAGS += -I$(MATLAB_ROOT)\simulink\include 
CFLAGS += -I$(MATLAB_ROOT)\rtw\c\src
CFLAGS += -I$(MATLAB_ROOT)\rtw\c\src\ext_mode\common
 
CFLAGS += -I$(MATLAB_ROOT)\rtw\c\src\ext_mode\custom 
CFLAGS += -I$(MATLAB_ROOT)\toolbox\coder\rtiostream\src\utils
CFLAGS += -I$(MATLAB_ROOT)\toolbox\coder\rtiostream\rtiostreamtcpip
CFLAGS += -I$(MATLAB_ROOT)\rtw\c\src\common 
CFLAGS += -I$(RTDACZYNQ_TOOLBOX_PATH)\blocks\c\ 

CFLAGS += -I$(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\include

CFLAGS += -DSAMPLING_FREQ=$(FREQUENCY) 
CFLAGS += -DUSE_RTMODEL -DMODEL=$(MODEL) -DRT -DNUMST=$(NUMST)  
CFLAGS += -DMODEL_STR="""$(MODEL)"""   
CFLAGS += -DTID01EQ=$(TID01EQ) -DNCSTATES=$(NCSTATES)  -DMT=0 -DHAVESTDIO 
CFLAGS += -DMAT_FILE=1  -DONESTEPFCN=1 -DTERMFCN=1  
CFLAGS += -DMULTI_INSTANCE_CODE=0  -DCLASSIC_INTERFACE=0 -DEXT_MODE -Wall -O0 -g3 -fmessage-length=0
CFLAGS += -c 
CFLAGS += -DRT_MODEL=RT_MODEL_$(MODEL)
CFLAGS += -DWAIT_FOR_START_PACKAGE=$(WAIT_FOR_START_PACKAGE)


MODULES_OBJ = $(MODULES:.c=.o)

OBJS  = $(MODEL).o 
OBJS += kk_rt_main.o  
OBJS += $(MODULES_OBJ)

OBJS += updown.o  
OBJS += ext_svr.o ext_work.o ext_svr_custom_transport.o unlink.o
OBJS += uZSimulink.o

# Libraries:
$(MODEL): $(OBJS) libs $(MAKEFILE_NAME)
	@echo ### Linking $(ELF_NAME)...
	$(LD) $(LDFLAGS) -o $(ELF_NAME) $(OBJS) -lm -Wl,--start-group,-lxil,-llocalfreertos,-lgcc,-lc,--end-group 
	@echo ### Created executable  $(ELF_NAME)

# Libraries:

INCLUDES  = -I$(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\local
INCLUDES += -I$(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\libsrc\freertos_zynq_v1_02_a\src
INCLUDES += -I$(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\include
LIBXIL = libxil.a
LOCALLIBFREERTOS = liblocalfreertos.a
COMPILER_FLAGS=  -O2 -c
EXTRA_COMPILER_FLAGS=-g -DUSE_AMP=1

FREERTOS_C_FILES := $(wildcard $(RTDACZYNQ_TOOLBOX_PATH)/AMPFreeRTOS_bsp/ps7_cortexa9_1/libsrc/freertos_zynq_v1_02_a/src/*.c)
FREERTOS_S_FILES := $(wildcard $(RTDACZYNQ_TOOLBOX_PATH)/AMPFreeRTOS_bsp/ps7_cortexa9_1/libsrc/freertos_zynq_v1_02_a/src/*.s)
FREERTOS_OBJ_FILES := $(FREERTOS_C_FILES:.c=.o) $(FREERTOS_S_FILES:.s=.o)

# Replace \ by / - required by the notdir function
OUTS1 = ${subst \,/,$(FREERTOS_OBJ_FILES)}
OUTS = $(notdir $(OUTS1))
FREERTOS_OBJ_LOCAL_FILES = $(notdir $(FREERTOS_OBJ_FILES))


libs: $(FREERTOS_OBJ_LOCAL_FILES) $(MAKEFILE_NAME)
	@echo ### libs...
	@echo "### FREERTOS_C_FILES:  $(FREERTOS_C_FILES)"
	@echo "### OUTS:       $(OUTS)"
	@echo "### TFINAL:     $(TFINAL)"
	@echo "### TPERIOD:    $(TPERIOD)"
	@echo "### FREQUENCY:  $(FREQUENCY)"
	@echo "### MAKEFILE_NAME:  $(MAKEFILE_NAME)"


	$(AR) -r ${LOCALLIBFREERTOS} ${OUTS}

%.o : $(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\libsrc\freertos_zynq_v1_02_a\src\%.c $(MAKEFILE_NAME)
	@echo ### Compiling FreeRTOS (*.c) $<
	$(CC) -DSAMPLING_FREQ=$(FREQUENCY) $(COMPILER_FLAGS) $(EXTRA_COMPILER_FLAGS) $(INCLUDES) $<	

%.o : $(RTDACZYNQ_TOOLBOX_PATH)\AMPFreeRTOS_bsp\ps7_cortexa9_1\libsrc\freertos_zynq_v1_02_a\src\%.s $(MAKEFILE_NAME)
	@echo ### Compiling FreeRTOS (*.s) $<
	$(CC) $(COMPILER_FLAGS) $(EXTRA_COMPILER_FLAGS) $(INCLUDES) $<	


%.o: $(RTDACZYNQ_TOOLBOX_PATH)\src\%.c $(MAKEFILE_NAME)
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o: $(MATLAB_ROOT)\work\%.c $(MAKEFILE_NAME)
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o: $(MATLAB_ROOT)\rtw\c\src\%.c $(MAKEFILE_NAME)
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o: $(MATLAB_ROOT)\rtw\c\src\ext_mode\common\%.c $(MAKEFILE_NAME)
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o: $(MATLAB_ROOT)\rtw\c\src\rtiostream\rtiostreamtcpip\%.c $(MAKEFILE_NAME)  
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o: $(RTDACZYNQ_TOOLBOX_PATH)\blocks\c\%.c $(MAKEFILE_NAME)  
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $< 

%.o : %.c $(MAKEFILE_NAME) 
	@echo ### Compiling $<
	$(CC) $(CFLAGS) $<

# Libraries:







