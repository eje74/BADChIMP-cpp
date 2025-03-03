# BADChIMP-cpp
## Build and compile 
**Linux:** Make sure that the openMPI libraries, including libopenmpi-dev, are installed. Run `./make.sh <name_of_folder_with_main_file>`, in the main folder, `std_case` is built if no argument to make.sh is given.  
This script will make a ```build``` folder, run ```cmake``` from that folder and then run ```make```. This can also be done by hand:
```shell
/BADChIMP-cpp$ mkdir <build-folder-name>
/BADChIMP-cpp$ cd <build-folder-name>
/BADChIMP-cpp$ cmake -DLBMAIN:STRING="<name_of_folder_with_main_file>" ./..
/BADChIMP-cpp$ make
``` 

**Windows:** Make sure that open [MPI is installed](https://docs.microsoft.com/en-us/archive/blogs/windowshpc/how-to-compile-and-run-a-simple-ms-mpi-program). Download and run `msmpisetup.exe` and `msmpisdk.msi`.  Install [cmake for Windows](https://cmake.org/). Run cmake from root directory to generate Visual Studio C++ project, or simply use VSCode.

