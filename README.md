# BADChIMP-cpp

---

## &#x1F411;&#x1F411;  Cloning the Repository (SSH, Linux CLI)  
This repository is cloned using SSH, which allows secure access without repeated login prompts.

---

### ✅ Step 1: Check for an existing SSH key
```bash
ls -la ~/.ssh
```  
Look for key pairs such as:
```bash
id_ed25519 and id_ed25519.pub (recommended)
id_rsa and id_rsa.pub (older)
```
If they exist, you can reuse them.

### ✅ Step 2: Create a new SSH key (if needed)
If no key exists, generate one:
```bash
ssh-keygen -t ed25519 -C "your.email@domain.com"
```
Press Enter to accept the default file location
Optionally set a passphrase for additional security

### ✅ Step 3: Add your SSH key to GitHub
Copy your public key:
```bash
cat ~/.ssh/id_ed25519.pub
```
This will printout the key to screen, mark it and copy it from there.  
Then in GitHub:
* Go to Profile Drop-down Menu →  Settings → SSH and GPG keys
* Click New SSH key
* Paste the key and save

### ✅ Step 4: Test the SSH connection
```bash
ssh -T git@github.com
```
The first time you connect, type yes to confirm the host fingerprint
A success message confirms authentication

### ✅ Step 5: Clone the repository
```bash
git clone git@github.com:eje74/BADChIMP-cpp.git
cd BADChIMP-cpp
```
---


## Build and compile 
**Linux:** Make sure that the openMPI libraries, including libopenmpi-dev, are installed. Run `./make.sh <name_of_folder_with_main_file>`, in the main directory, `std_case` is built if no argument to make.sh is given.  
This script will make a ```build``` folder, run ```cmake``` from that folder and then run ```make```. This can also be done by hand:
```shell
/BADChIMP-cpp$ mkdir <build-folder-name>
/BADChIMP-cpp$ cd <build-folder-name>
/BADChIMP-cpp$ cmake -DLBMAIN:STRING="<name_of_folder_with_main_file>" ./..
/BADChIMP-cpp$ make
``` 

**Windows:** Make sure that open [MPI is installed](https://docs.microsoft.com/en-us/archive/blogs/windowshpc/how-to-compile-and-run-a-simple-ms-mpi-program). Download and run `msmpisetup.exe` and `msmpisdk.msi`.  Install [cmake for Windows](https://cmake.org/). Run cmake from main directory to generate Visual Studio C++ project, or simply use VSCode.

