# Makefile

# 编译器
NVCC = nvcc
GXX = g++

# 编译选项
NVCCFLAGS = -g -O2
GXXFLAGS = -g -O2 -std=c++17

# 源文件
CU_SRCS = main-montgomery-GPU.cu
CPP_SRCS = montgomery.cpp

# 目标文件
CU_OBJS = $(CU_SRCS:.cu=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)

TARGET = main-montgomery-GPU

.PHONY: all clean

all: $(TARGET)

# 编译 CUDA 文件，只编译，不链接
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# 编译 C++ 文件，只编译，不链接
%.o: %.cpp
	$(GXX) $(GXXFLAGS) -c $< -o $@

# 链接所有目标文件，使用 g++ 链接保证支持 C++ 标准库和特性
$(TARGET): $(CU_OBJS) $(CPP_OBJS)
	$(GXX) $(GXXFLAGS) $^ -o $@ -L/usr/local/cuda/lib64 -lcudart

clean:
	rm -f $(CU_OBJS) $(CPP_OBJS) $(TARGET)
