Run
```
mkdir -p bin build doc include lib snippets src test
```
to make the file structure

* **bin**: The output executables go here, both for the app and for any tests and
snippets.
* **build**: This folder contains all object files, and is removed on a `clean`.
* **config**: development and production config files in here
* **doc**: Any notes, like my assembly notes and configuration files, are here.
* **include**: All project header files. All necessary third-party header files that
do not exist under `/usr/local/include` are also placed here.
* **lib**: Any libs that get compiled by the project, third party or any needed in
development. Prior to deployment, third party libraries get moved to
`/usr/local/lib` where they belong, leaving the project clean enough to compile
on our Linux deployment servers. I really use this to test different library
versions than the standard.
* **snippets**: I often write smaller classes or files to test technologies or
ideas, and keep them around for future reference. They go here, where they do
not dilute the real application’s files, but can still be found later.
* **src**: The application and only the application’s source files.
* **test**: All test code files. You do write tests, no?  

This is taken from [Hiltmon's blog](https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/)
