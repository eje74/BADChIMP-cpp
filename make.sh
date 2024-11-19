# pass as argument the name of the directory where the main file is located
# std_case will be used if directory not found. e.g.
# make twophase
# make turbine
# etc.

# Use a common ../make.sh file or do customized build here.
function system {
  "$@"
  if [ $? -ne 0 ]; then
    echo "make.sh: unsuccessful command $@"
    echo "abort!"
    exit 1
  fi
}
dest=build
if [ ! -d $dest ]; then
  mkdir $dest
fi

system cd $dest
system cmake -DLBMAIN:STRING="$@" ../
system make
echo "executable in: build/bin/bdchmp"
echo "to build other cases run: e.g. make <folder-name-where-main-file-is-located> (e.g. twophase)"
